#ifndef FWDPP_TS_MARGINAL_TREE_FUNCTIONS_CHILDREN_HPP
#define FWDPP_TS_MARGINAL_TREE_FUNCTIONS_CHILDREN_HPP

#include "../marginal_tree.hpp"

namespace fwdpp
{
    namespace ts
    {
        class child_iterator
        {
          private:
            bool left_to_right;
            TS_NODE_INT current_child;
            std::vector<TS_NODE_INT>::const_iterator sib_begin, sib_end;

            inline TS_NODE_INT
            init_current_child(const marginal_tree& m, TS_NODE_INT n)
            {
                if (n == TS_NULL_NODE)
                    {
                        throw std::invalid_argument("node is NULL");
                    }
                if (static_cast<std::size_t>(n) >= m.left_child.size())
                    {
                        throw std::invalid_argument("node id is out of range");
                    }
                if (left_to_right)
                    {
                        return m.left_child[n];
                    }
                return m.right_child[n];
            }

          public:
            child_iterator(const marginal_tree& m, TS_NODE_INT n,
                           bool left_to_right_traversal)
                : left_to_right(left_to_right_traversal),
                  current_child(init_current_child(m, n)),
                  sib_begin(left_to_right ? begin(m.right_sib)
                                          : begin(m.left_sib)),
                  sib_end(left_to_right ? end(m.right_sib) : end(m.left_sib))
            {
            }

            inline TS_NODE_INT
            operator()()
            {
                auto c = current_child;
                if (sib_begin + current_child >= sib_end)
                    {
                        throw std::runtime_error("child iteration error");
                    }
                current_child = *(sib_begin + current_child);
                return c;
            }

            template <typename F>
            inline bool
            operator()(const F& f)
            {
                auto c = this->operator()();
                bool rv = (c != TS_NULL_NODE);
                if (rv)
                    {
                        f(c);
                    }
                return rv;
            }
        };

        template <typename F>
        inline void
        process_children(const marginal_tree& m, TS_NODE_INT n,
                         bool left_to_right_traversal, const F& f)
        {
            child_iterator ci(m, n, left_to_right_traversal);
            while (ci(f))
                {
                }
        }

        inline std::vector<TS_NODE_INT>
        get_children(const marginal_tree& m, TS_NODE_INT n,
                     bool left_to_right_traversal)
        {
            std::vector<TS_NODE_INT> rv;
            process_children(m, n, left_to_right_traversal,
                             [&rv](const TS_NODE_INT c) { rv.push_back(c); });
            return rv;
        }

        inline int
        num_children(const marginal_tree& m, TS_NODE_INT n)
        {
            int nc = 0;
            process_children(m, n, true,
                             [&nc](const TS_NODE_INT /*c*/) { ++nc; });
            return nc;
        }
    } // namespace ts
} // namespace fwdpp

#endif
