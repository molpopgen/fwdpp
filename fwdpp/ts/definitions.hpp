#ifndef FWDPP_TS_DEFINITIONS_HPP
#define FWDPP_TS_DEFINITIONS_HPP

#include <cstdint>

/// \namespace fwdpp::ts Tree sequence \cite Kelleher2018-fu support

namespace fwdpp
{
    namespace ts
    {
        /// Integer type for table indexes
        using table_index_t [[deprecated("Use member typedefs instead")]] = std::int32_t;
        /// NULL index value
        constexpr std::int32_t NULL_INDEX
            [[deprecated("Use member constexpr value null instead")]]
            = -1;
        /// Convention for the ancestral state of a site
        constexpr std::int8_t default_ancestral_state = 0;
        /// Convention for the derived state of a site
        constexpr std::int8_t default_derived_state = 1;
    } // namespace ts
} // namespace fwdpp

#endif
