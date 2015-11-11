#ifndef LIBGRAPHLET_SIMILARITYGPU_HPP
#define LIBGRAPHLET_SIMILARITYGPU_HPP

#include <vector>
#include <boost/compute/core.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <orca/Orca.hpp>

namespace libgraphlet {
	void similarityGPU(
		const orca::Orca &oa,
		const orca::Orca &ob,
		boost::numeric::ublas::matrix<float> &sim,
		const boost::compute::device &device = boost::compute::system::default_device()
	);
}

#endif
