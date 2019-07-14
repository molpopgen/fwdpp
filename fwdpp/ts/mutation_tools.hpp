#ifndef FWDPP_TS_MUTATION_TOOLS_HPP__
#define FWDPP_TS_MUTATION_TOOLS_HPP__

#include "site.hpp"

namespace fwdpp
{
    namespace ts
    {
        struct new_variant_record
        {
            site s;
            std::size_t key;
            bool neutral;
            std::int8_t derived_state;
            new_variant_record(double position, std::int8_t ancestral_state,
                               std::size_t index, bool n, std::int8_t d)
                : s{ position, ancestral_state }, key{ index }, neutral{ n },
                  derived_state{ d }
            /// \version 0.8.0 Added to library
            {
            }
        };
    } // namespace ts
} // namespace fwdpp

#endif
