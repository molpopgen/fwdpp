#ifndef FWDPP_TS_MUTATION_RECORD_HPP
#define FWDPP_TS_MUTATION_RECORD_HPP

#include <cstdint>
#include "types/mutation_record.hpp"

namespace fwdpp
{
    namespace ts
    {
        /// 32-bit mutation_record
        using mutation_record = types::mutation_record<std::int32_t>;
    } // namespace ts
} // namespace fwdpp

#endif
