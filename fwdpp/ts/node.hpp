// Copied from fwdpy11_arg_example.
// Author: KRT
// License: GPL3+
#ifndef FWDPP_TS_NODE_HPP
#define FWDPP_TS_NODE_HPP

#include <cstdint>

namespace fwdpp
{
    namespace ts
    {
        struct node
        {
            std::int32_t population;
            double generation;
        };
    }
}
#endif
