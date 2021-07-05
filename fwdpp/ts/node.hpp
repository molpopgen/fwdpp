// Copied from fwdpy11_arg_example.
// Author: KRT
// License: GPL3+
#ifndef FWDPP_TS_NODE_HPP
#define FWDPP_TS_NODE_HPP

#include <cstdint>
#include "types/node.hpp"

namespace fwdpp
{
    namespace ts
    {
        /// 32-bit node
        using node = types::node<std::int32_t>;
    } // namespace ts
} // namespace fwdpp
#endif
