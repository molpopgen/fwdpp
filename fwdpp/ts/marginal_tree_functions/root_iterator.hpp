#ifndef FWDPP_TS_MARGINAL_TREE_ROOT_ITERATOR_HPP
#define FWDPP_TS_MARGINAL_TREE_ROOT_ITERATOR_HPP

#include <stdexcept>
#include "../marginal_tree.hpp"

namespace fwdpp
{
    namespace ts
    {
        class root_iterator
        /// Allows iteration over tree roots
        /// \note Tied to lifetime of a marginal_tree.
        {
          private:
            TS_NODE_INT current_root;
            std::vector<TS_NODE_INT>::const_iterator rsib_beg, rsib_end;

          public:
            explicit root_iterator(const marginal_tree& m)
                : current_root(m.left_root), rsib_beg(begin(m.right_sib)),
                  rsib_end(end(m.right_sib))
            {
                if (current_root == TS_NULL_NODE)
                    {
                        throw std::invalid_argument("root list is empty");
                    }
                if (rsib_beg == rsib_end)
                    {
                        throw std::invalid_argument(
                            "empty list of right sibs");
                    }
            }

            inline TS_NODE_INT
            operator()()
            /// TS_NULL_NODE signals end of iteration
            {
                auto croot = current_root;
                if (rsib_beg + croot >= rsib_end)
                    {
                        throw std::runtime_error("root iteration error");
                    }
                current_root = *(rsib_beg + croot);
                return croot;
            }
        };
    } // namespace ts
} // namespace fwdpp
#endif
