#ifndef FWDPP_TS_SITE_HPP
#define FWDPP_TS_SITE_HPP

#include <tuple>
#include <cstdint>

namespace fwdpp
{
    namespace ts
    {
        struct site
        /// Entry in a site table
        /// \version 0.8.0 Added to llibrary
        {
            double position;
            std::int8_t ancestral_state;
        };

        inline bool
        operator<(const site& a, const site& b)
        {
            return a.position < b.position;
        }

        inline bool
        operator==(const site& a, const site& b)
        {
            return std::tie(a.position, a.ancestral_state)
                   == std::tie(b.position, b.ancestral_state);
        }
    } // namespace ts
} // namespace fwdpp

#endif

