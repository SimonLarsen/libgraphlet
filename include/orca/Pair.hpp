#ifndef ORCA_PAIR_HPP
#define ORCA_PAIR_HPP

#include <algorithm>

namespace orca {
	class Pair {
		public:
			Pair() : a(0), b(0) { }

			Pair(int a0, int b0) {
				a = std::min(a0,b0);
				b = std::max(a0,b0);
			}

			int a, b;
	};

	inline bool operator<(const Pair &x, const Pair &y) {
		if (x.a==y.a) return x.b<y.b;
		else return x.a<y.a;
	}

	inline bool operator==(const Pair &x, const Pair &y) {
		return x.a==y.a && x.b==y.b;
	}

	struct HashPair {
		size_t operator()(const Pair &x) const {
			return (x.a<<8) ^ (x.b<<0);
		}
	};
}

#endif
