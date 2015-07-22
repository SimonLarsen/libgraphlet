#ifndef ORCA_SIMILARITY_HPP
#define ORCA_SIMILARITY_HPP

#include <vector>
#include <boost/compute/core.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <orca/Orca.hpp>

namespace orca {
	void similarity(
		const Orca &oa,
		const Orca &ob,
		boost::numeric::ublas::matrix<float> &sim,
		const boost::compute::device &device = boost::compute::system::default_device()
	);
}

#endif
