#include <orca/Similarity.hpp>

#include <iostream>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <thread>
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
		if(orbits > 1) {
			for(size_t k = 0; k < orbits; ++k) {
				weights[k] = 1.0f - log(AFFECTED[k]) / log(orbits);
			}
		} else {
			weights[0] = 1.0f;
		}
		float weights_sum = std::accumulate(weights.begin(), weights.end(), 0.0f);

		sim.resize(na, nb);

		auto start_time = std::chrono::system_clock::now();

		size_t nthreads = std::thread::hardware_concurrency();
		std::vector<std::thread> threads;

		for(size_t p = 0; p < nthreads; ++p) {
			threads.emplace_back(std::thread(
				[&, p](){
					for(size_t i = p; i < na; i += nthreads) {
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
			));
		}

		for(auto &t : threads) {
			t.join();
		}

		auto end_time = std::chrono::system_clock::now();

		std::cerr << "time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time).count() << " ms" << std::endl;
	}
}
