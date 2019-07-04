#ifndef FWDPP_TS_MARGINAL_TREE_NUM_ROOTS_HPP
#define FWDPP_TS_MARGINAL_TREE_NUM_ROOTS_HPP

#include <stdexcept>
#include "../marginal_tree.hpp"

namespace fwdpp
{
    namespace ts
    {
        inline int
        num_roots(const marginal_tree& m)
        /// Return number of roots
        {
            auto lr = m.left_root;
            if (lr == TS_NULL_NODE)
                {
                    throw std::invalid_argument("left_root is NULL");
                }
            int nroots = 0;
            while (lr != TS_NULL_NODE)
                {
                    ++nroots;
                    lr = m.right_sib[lr];
                }
            return nroots;
        }
    } // namespace ts
} // namespace fwdpp

#endif
