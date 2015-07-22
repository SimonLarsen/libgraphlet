#ifndef ORCA_TRIPLE_HPP
#define ORCA_TRIPLE_HPP

#include <algorithm>

namespace orca {
	class Triple {
		public:
			int a, b, c;

			Triple(int a0, int b0, int c0) {
				a = a0;
				b = b0;
				c = c0;
				if (a > b) std::swap(a,b);
				if (b > c) std::swap(b,c);
				if (a > b) std::swap(a,b);
			}
	};

	inline bool operator<(const Triple &x, const Triple &y) {
		if (x.a==y.a) {
			if (x.b==y.b) return x.c<y.c;
			else return x.b<y.b;
		} else return x.a<y.a;
	}

	inline bool operator==(const Triple &x, const Triple &y) {
		return x.a==y.a && x.b==y.b && x.c==y.c;
	}

	struct HashTriple {
		size_t operator()(const Triple &x) const {
			return (x.a<<16) ^ (x.b<<8) ^ (x.c<<0);
		}
	};
}

#endif
