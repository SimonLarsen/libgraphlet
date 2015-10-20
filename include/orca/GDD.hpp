#ifndef ORCA_GDD_HPP
#define ORCA_GDD_HPP 

#include <vector>
#include <map>
#include <orca/Orca.hpp>

namespace orca {
	typedef std::vector<std::map<size_t, float>> GDD;

	void gdd(
		const Orca &orca,
		GDD &gdd,
		bool normalize = true
	);

	void gdd_agreement(
		const GDD &a,
		const GDD &b,
		std::vector<float> &out
	);
}

#endif
