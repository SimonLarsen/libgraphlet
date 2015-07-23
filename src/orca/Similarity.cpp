#include <orca/Similarity.hpp>

#include <iostream>
#include <algorithm>
#include <boost/log/trivial.hpp>
#include <orca/OrcaException.hpp>

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
	void similarity(
		const Orca &oa,
		const Orca &ob,
		boost::numeric::ublas::matrix<float> &sim
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

		sim.resize(na, nb);

		for(size_t i = 0; i < na; ++i) {
			for(size_t j = 0; j < nb; ++j) {
				float D = 0.0f;
				for(size_t k = 0; k < orbits; ++k) {
					int64_t aik = oa.getOrbits()(i, k);
					int64_t bjk = ob.getOrbits()(j, k);

					float num = fabs(log(aik + 1.0f) - log(bjk + 1.0f));
					float denom = log(std::max(aik, bjk) + 2.0f);
					D += weights[k] * num / denom;
				}

				sim(i, j) = 1.0f - D / weights_sum;
			}
		}
	}
}
