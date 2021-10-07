#pragma once

#include <tuple>
#include <cstdint>
#include <limits>

namespace fwdpp
{
    namespace ts
    {
        namespace types
        {
            struct site
            /// Entry in a site table
            /// \version 0.8.0 Added to llibrary
            {
                double position;
                std::int8_t ancestral_state;

                site()
                    : position{std::numeric_limits<double>::quiet_NaN()},
                      ancestral_state{-1}
                {
                }

                site(double position, std::int8_t ancestral_state)
                    : position{position}, ancestral_state{ancestral_state}
                {
                }
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
        }
    }
}
