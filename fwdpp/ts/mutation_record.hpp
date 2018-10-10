#ifndef FWDPP_TS_MUTATION_RECORD_HPP
#define FWDPP_TS_MUTATION_RECORD_HPP

#include <cstdint>

namespace fwdpp
{
	namespace ts
	{
		struct mutation_record
		{
            /// The node to which the mutation is
            /// currently simplified
			std::int32_t node;
            /// The index of the mutation in the
            /// population's mutation container
			std::size_t key;
		};
	}
}

#endif
