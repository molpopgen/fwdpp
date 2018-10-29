// Copied from fwdpy11_arg_example.
// Author: KRT
// License: GPL3+
#ifndef FWDPP_TS_NODE_HPP
#define FWDPP_TS_NODE_HPP

#include <tuple>
#include <cstdint>

namespace fwdpp
{
    namespace ts
    {
        struct node
        /// A node in a tree sequence
        ///  \version 0.7.0 Added to fwdpp
        {
            /// Location of the node.
            /// Used for models of discrete population structure
            std::int32_t population;
            /// Birth time of the node.
            double time;
        };

        bool
        operator==(const node& a, const node& b)
        {
            return std::tie(a.time, a.population)
                   == std::tie(b.time, b.population);
        }
    } // namespace ts
} // namespace fwdpp
#endif
