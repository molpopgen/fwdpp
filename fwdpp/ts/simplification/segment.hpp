#ifndef FWDPP_TS_SIMPLIFICATION_SEGMENT_HPP__
#define FWDPP_TS_SIMPLIFICATION_SEGMENT_HPP__

#include <fwdpp/ts/definitions.hpp>
#include <stdexcept>

namespace fwdpp
{
    namespace ts
    {
        struct segment
        {
            double left, right;
            TS_NODE_INT node;
            segment(double l, double r, TS_NODE_INT n)
                : left{ l }, right{ r }, node{ n }
            {
                if (right <= left)
                    {
                        throw std::invalid_argument("right must be > left");
                    }
            }
        };

    } // namespace ts
} // namespace fwdpp

#endif
