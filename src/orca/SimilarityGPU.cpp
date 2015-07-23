#include <orca/SimilarityGPU.hpp>

#include <iostream>
#include <algorithm>
#include <boost/compute/container/mapped_view.hpp>
#include <boost/compute/container/vector.hpp>
#include <orca/OrcaException.hpp>
#include "Kernels.hpp"

namespace compute = boost::compute;

namespace {
	const int AFFECTED[73] = {
		1, 2, 2, 2, 3, 4, 3, 3, 4, 3,
		4, 4, 4, 4, 3, 4, 6, 5, 4, 5,
		6, 6, 4, 4, 4, 5, 7, 4, 6, 6,
		7, 4, 6, 6, 6, 5, 6, 7, 7, 5,
		7, 6, 7, 6, 5, 5, 6, 8, 7, 6,
		6, 8, 6, 9, 5, 6, 4, 6, 6, 7,
		8, 6, 6, 8, 7, 6, 7, 7, 8, 5,
		6, 6, 4
	};
}

namespace orca {
	void similarityGPU(
		const Orca &oa,
		const Orca &ob,
		boost::numeric::ublas::matrix<float> &sim,
		const compute::device &device
	) {
		if(oa.graphletSize() != ob.graphletSize()) {
			throw OrcaException("Orca instances not of same size.");
		}
		int graphlet_size = oa.graphletSize();
		const size_t orbits = ORBITS[graphlet_size];
		const size_t na = oa.getOrbits().size1();
		const size_t nb = ob.getOrbits().size1();

		std::vector<float> weights(orbits);
		for(size_t k = 0; k < orbits; ++k) {
			weights[k] = 1.0f - log(AFFECTED[k]) / log(orbits);
		}
		float weights_sum = std::accumulate(weights.begin(), weights.end(), 0.0f);

		compute::context context(device);
		compute::command_queue queue(context, device);

		compute::program program = compute::program::create_with_source(
			kernel_orca_similarity_source,
			context
		);

		try {
			program.build();
		} catch(compute::opencl_error &e) {
			std::cerr << program.build_log();
			throw;
		}

		compute::kernel kernel = program.create_kernel("orca_compute_similarity");

		compute::vector<cl_uint> buf_weights(orbits, context);
		compute::vector<cl_float> buf_sim(na * nb, context);

		compute::mapped_view<cl_long> map_a(&(oa.getOrbits().data()[0]), na*orbits, context);
		compute::mapped_view<cl_long> map_b(&(ob.getOrbits().data()[0]), nb*orbits, context);

		int arg = 0;
		kernel.set_arg(arg++, (cl_uint)na);
		kernel.set_arg(arg++, (cl_uint)nb);
		kernel.set_arg(arg++, (cl_uint)orbits);
		kernel.set_arg(arg++, (cl_uint)weights_sum);
		kernel.set_arg(arg++, buf_weights);
		kernel.set_arg(arg++, map_a.get_buffer());
		kernel.set_arg(arg++, map_b.get_buffer());
		kernel.set_arg(arg++, buf_sim);

		compute::copy(
			weights.begin(),
			weights.end(),
			buf_weights.begin(),
			queue
		);
		queue.enqueue_1d_range_kernel(kernel, 0, 1024, 0);

		sim.resize(na, nb);
		compute::copy(
			buf_sim.begin(),
			buf_sim.end(),
			&(sim.data()[0]),
			queue
		);

		queue.finish();
	}
}
