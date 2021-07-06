#ifndef FWDPP_TS_MARGINAL_TREE_FUNCTIONS_CHILDREN_HPP
#define FWDPP_TS_MARGINAL_TREE_FUNCTIONS_CHILDREN_HPP

#include "../marginal_tree.hpp"
#include "../types/generate_null_id.hpp"

namespace fwdpp
{
    namespace ts
    {
        template <typename SignedInteger> class child_iterator
        /// \brief Faciliate traversal of a node's children.
        /// \headerfile fwdpp/ts/marginal_tree<SignedInteger>_functions/children.hpp
        {
          private:
            bool left_to_right;
            SignedInteger current_child;
            typename std::vector<SignedInteger>::const_iterator sib_begin, sib_end;

            inline SignedInteger
            init_current_child(const marginal_tree<SignedInteger>& m, SignedInteger n)
            {
                if (n == types::generate_null_id<SignedInteger>())
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
            child_iterator(const marginal_tree<SignedInteger>& m, SignedInteger n,
                           bool left_to_right_traversal)
                : left_to_right(left_to_right_traversal),
                  current_child(init_current_child(m, n)),
                  sib_begin(left_to_right ? begin(m.right_sib) : begin(m.left_sib)),
                  sib_end(left_to_right ? end(m.right_sib) : end(m.left_sib))
            /// \param m A fwdpp::ts::marginal_tree<SignedInteger>
            /// \param n Index of a fwdpp::ts::node
            /// \param left_to_right_traversal.  If `true`, traverse children left-to-right. Else, the opposite.
            {
            }

            inline SignedInteger
            operator()()
            /// Advance to the next child.
            ///
            /// Returns fwdpp::ts::types::generate_null_id<SignedInteger>() when no more children remain.
            {
                static_assert(types::generate_null_id<SignedInteger>() < 0,
                              "types::generate_null_id<SignedInteger>() < 0 is false, "
                              "so something needs to change here");
                auto c = current_child;
                if (c < 0)
                    {
                        return c;
                    }
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
            /// Apply a function to each child
            /// \param f A function equivalent to void (*process_child)(SignedInteger)
            ///
            /// Returns false to signify end of iteration.
            {
                auto c = this->operator()();
                bool rv = (c != types::generate_null_id<SignedInteger>());
                if (rv)
                    {
                        f(c);
                    }
                return rv;
            }
        };

        template <typename SignedInteger, typename F>
        inline void
        process_children(const marginal_tree<SignedInteger>& m, SignedInteger n,
                         bool left_to_right_traversal, const F& f)
        /// \brief Apply a function to children of node \a n
        /// \param m A fwdpp::ts::marginal_tree<SignedInteger>
        /// \param n Index a fwdpp::ts::node
        /// \param left_to_right_traversal.  If `true`, traverse children left-to-right. Else, the opposite.
        /// \param f A function equivalent to void (*process_child)(SignedInteger)
        {
            child_iterator<SignedInteger> ci(m, n, left_to_right_traversal);
            while (ci(f))
                {
                }
        }

        template <typename SignedInteger>
        inline std::vector<SignedInteger>
        get_children(const marginal_tree<SignedInteger>& m, SignedInteger n,
                     bool left_to_right_traversal)
        /// \brief Get all children of a node
        /// \param m A fwdpp::ts::marginal_tree<SignedInteger>
        /// \param n Index a fwdpp::ts::node
        /// \param left_to_right_traversal.  If `true`, traverse children left-to-right. Else, the opposite.
        /// \return std::vector<fwdpp::ts::SignedInteger>
        {
            std::vector<SignedInteger> rv;
            process_children(m, n, left_to_right_traversal,
                             [&rv](const SignedInteger c) { rv.push_back(c); });
            return rv;
        }

        template <typename SignedInteger>
        inline int
        num_children(const marginal_tree<SignedInteger>& m, SignedInteger n)
        /// \brief Return number of children of a node
        /// \param m A fwdpp::ts::marginal_tree<SignedInteger>
        /// \param n Index a fwdpp::ts::node
        /// \return int
        {
            int nc = 0;
            process_children(m, n, true, [&nc](const SignedInteger /*c*/) { ++nc; });
            return nc;
        }
    } // namespace ts
} // namespace fwdpp

#endif
