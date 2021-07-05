//Copied from fwdpy11_arg_example.
//Author: KRT
//License: GPL3+
#ifndef FWDPP_TS_EDGE_HPP
#define FWDPP_TS_EDGE_HPP

#include <cstdint>
#include "types/edge.hpp"

namespace fwdpp
{
    namespace ts
    {
        /// 32-bit edge
        using edge = types::edge<std::int32_t>;
    } // namespace ts
} // namespace fwdpp

#endif
