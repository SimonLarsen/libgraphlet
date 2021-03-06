#include <iostream>
#include <stdexcept>
#include <libgraphlet/GDD.hpp>

namespace libgraphlet {
	void gdd(
		const orca::Orca &orca,
		GDD &gdd,
		bool normalize
	) {
		int orbits = orca::ORBITS[orca.graphletSize()];

		gdd.clear();
		gdd.resize(orbits);

		auto &sig = orca.getOrbits();

		// Calculate degree distributions
		for(auto it1 = sig.begin1(); it1 != sig.end1(); ++it1) {
			size_t j = 0;
			for(auto it = it1.begin(); it != it1.end(); ++it, ++j) {
				auto k = *it;

				if(k == 0) continue;

				if(gdd[j].find(k) == gdd[j].end()) {
					gdd[j][k] = 0.0f;
				}

				gdd[j][k] += 1.0f;
			}
		}

		if(normalize) {
			// Scale by k
			for(auto &v : gdd) {
				for(auto &it : v) {
					if(it.first > 0) {
						it.second /= (float)it.first;
					}
				}
			}

			// Normalize to [0,1]
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

	/**
	 * Calculates the GDD-agreement vector of two GDDs.
	 * Note: GDDs must be normalized!
	 */
	void gdd_agreement(
		const GDD &a,
		const GDD &b,
		std::vector<float> &out
	) {
		if(a.size() != b.size()) {
			throw std::invalid_argument("GDDs not of same size");
		}

		size_t n = a.size();
		out.resize(n);

		for(size_t i = 0; i < n; ++i) {
			size_t a_max = (a[i].size() > 0 ? a[i].rbegin()->first : 0);
			size_t b_max = (b[i].size() > 0 ? b[i].rbegin()->first : 0);
			size_t m = std::max(a_max, b_max);

			float sum = 0.0f;
			for(size_t j = 1; j <= m; ++j) {
				float a_value = 0.0f;
				float b_value = 0.0f;

				if(a[i].find(j) != a[i].end()) {
					a_value = a[i].at(j);
				}
				if(b[i].find(j) != b[i].end()) {
					b_value = b[i].at(j);
				}

				sum += std::pow(a_value - b_value, 2.0f);
			}

			out[i] = 1.0f - std::sqrt(sum) / std::sqrt(2.0f);
		}
	}
}
