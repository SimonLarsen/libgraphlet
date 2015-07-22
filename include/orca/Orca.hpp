#ifndef ORCA_ORCA_HPP
#define ORCA_ORCA_HPP

#include <unordered_map>
#include <utility>
#include <boost/numeric/ublas/matrix.hpp>
#include <orca/Pair.hpp>
#include <orca/Triple.hpp>

namespace orca {
	typedef boost::numeric::ublas::matrix<int64_t> Signature;

	unsigned const int ORBITS[6] = { 0, 0, 1, 4, 15, 73 };

	class Orca {
		public:
			Orca(
				size_t n,
				const std::vector<std::pair<size_t,size_t>> &in_edges,
				unsigned int graphlet_size
			);
			void compute();
			const Signature &getOrbits() const;
			int graphletSize() const;

		private:
			void count4();
			void count5();

			bool adjacent(int x, int y) const;

			int common3_get(int a, int b, int c) const;
			int common2_get(int a, int c) const;

			int n, m;
			unsigned int graphlet_size;
			std::vector<int> deg;
			std::vector<Pair> edges;
			std::vector<std::vector<int>> adj;
			std::vector<std::vector<std::pair<int,int>>> inc;
			Signature orbit;

			std::unordered_map<Pair, int, HashPair> common2;
			std::unordered_map<Triple, int, HashTriple> common3;
	};
}

#endif
