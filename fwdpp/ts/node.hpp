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
        ///  \version 0.8.0 Rename ::population to ::deme
        {
            /// Location of the node.
            /// Used for models of discrete population structure
            std::int32_t deme;
            /// Birth time of the node.
            double time;
        };

        inline bool
        operator==(const node& a, const node& b)
        {
            return std::tie(a.time, a.deme) == std::tie(b.time, b.deme);
        }
    } // namespace ts
} // namespace fwdpp
#endif
