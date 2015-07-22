#ifndef ORCA_ORCAEXCEPTION_HPP
#define ORCA_ORCAEXCEPTION_HPP

#include <stdexcept>

namespace orca {
	class OrcaException : public std::runtime_error {
		public:
			OrcaException(const std::string &message ) : std::runtime_error(message) { }
	};
}

#endif
