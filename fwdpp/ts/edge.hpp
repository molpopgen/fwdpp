//Copied from fwdpy11_arg_example.
//Author: KRT
//License: GPL3+
#ifndef FWDPP_TS_EDGE_HPP
#define FWDPP_TS_EDGE_HPP

#include "definitions.hpp"

namespace fwdpp
{
    namespace ts
    {
        struct edge
        {
            double left, right;
            TS_NODE_INT parent, child;
        };
    }
}

#endif
