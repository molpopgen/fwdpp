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
        /// \brief An edge in a tree sequence
		///
        /// Edges define a transmission event
        /// of the genomic interval [left,right)
        /// from parent to child.
        ///
        ///  \version 0.7.0 Added to fwdpp
        {
            double left, right;
            TS_NODE_INT parent, child;
        };
    } // namespace ts
} // namespace fwdpp

#endif
