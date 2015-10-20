#include <iostream>
#include <orca/GDD.hpp>

namespace orca {
	void gdd(
		const Orca &orca,
		std::vector<std::map<size_t,float>> &gdd,
		bool normalize
	) {
		int orbits = ORBITS[orca.graphletSize()];

		gdd.clear();
		gdd.resize(orbits);

		auto &sig = orca.getOrbits();

		// Calculate djG(k)
		for(auto it1 = sig.begin1(); it1 != sig.end1(); ++it1) {
			size_t j = 0;
			for(auto it = it1.begin(); it != it1.end(); ++it, ++j) {
				auto k = *it;

				if(gdd[j].find(k) == gdd[j].end()) {
					gdd[j][k] = 0.0f;
				}

				gdd[j][k] += 1.0f;
			}
		}

		// Scale by k
		if(normalize) {
			for(auto &v : gdd) {
				for(auto &it : v) {
					if(it.first > 0) {
						it.second /= (float)it.first;
					}
				}
			}

			// Normalize
			for(auto &v : gdd) {
				float sum = 0.0f;
				for(auto &it : v) {
					sum += it.second;
				}

				for(auto &it : v) {
					it.second /= sum;
				}
			}
		}
	}
}
