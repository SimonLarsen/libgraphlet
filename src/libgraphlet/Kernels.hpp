#ifndef LIBGRAPHLET_KERNELS_HPP
#define LIBGRAPHLET_KERNELS_HPP

#include <boost/compute/utility/source.hpp>

namespace libgraphlet {
const char *kernel_orca_similarity_source = R"source(

__kernel void orca_compute_similarity(
	const uint na,
	const uint nb,
	const uint orbits,
	const float weights_sum,
	__constant float *weights,
	__global long *a,
	__global long *b,
	__global float *sim
) {
	size_t global_id = get_global_id(0);
	size_t global_size = get_global_size(0);

	for(size_t i = global_id; i < na; i += global_size) {
		for(size_t j = 0; j < nb; ++j) {
			// Calculate distance
			float D = 0.0f;
			for(int k = 0; k < orbits; ++k) {
				ulong aik = a[i*orbits + k];
				ulong bjk = b[j*orbits + k];

				float num = fabs(log(aik + 1.0f) - log(bjk + 1.0f));
				float denom = log(max(aik, bjk) + 2.0f);
				D += weights[k] * num / denom;
			}

			// Map distance to similarity
			sim[i*nb + j] = 1.0f - D / weights_sum;
		}
	}
}

)source";

}

#endif
