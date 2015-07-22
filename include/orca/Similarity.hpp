#ifndef ORCA_SIMILARITY_HPP
#define ORCA_SIMILARITY_HPP

#include <vector>
#include <orca/Orca.hpp>

namespace orca {
	void similarity(
		const Orca &oa,
		const Orca &ob,
		boost::numeric::ublas::matrix<float> &sim
	);
}

#endif
