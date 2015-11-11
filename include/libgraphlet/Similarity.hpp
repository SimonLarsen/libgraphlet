#ifndef LIBGRAPHLET_SIMILARITY_HPP
#define LIBGRAPHLET_SIMILARITY_HPP

#include <vector>
#include <orca/Orca.hpp>

namespace libgraphlet {
	void similarity(
		const orca::Orca &oa,
		const orca::Orca &ob,
		boost::numeric::ublas::matrix<float> &sim
	);
}

#endif
