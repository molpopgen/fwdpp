#ifndef FWDPP_TS_MARGINAL_TREE_FUNCTIONS_CHILDREN_HPP
#define FWDPP_TS_MARGINAL_TREE_FUNCTIONS_CHILDREN_HPP

#include "../marginal_tree.hpp"

namespace fwdpp
{
    namespace ts
    {
        class child_iterator
        /// \brief Faciliate traversal of a node's children.
        /// \headerfile fwdpp/ts/marginal_tree_functions/children.hpp
        {
          private:
            bool left_to_right;
            table_index_t current_child;
            std::vector<table_index_t>::const_iterator sib_begin, sib_end;

            inline table_index_t
            init_current_child(const marginal_tree& m, table_index_t n)
            {
                if (n == NULL_INDEX)
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
            child_iterator(const marginal_tree& m, table_index_t n,
                           bool left_to_right_traversal)
                : left_to_right(left_to_right_traversal),
                  current_child(init_current_child(m, n)),
                  sib_begin(left_to_right ? begin(m.right_sib)
                                          : begin(m.left_sib)),
                  sib_end(left_to_right ? end(m.right_sib) : end(m.left_sib))
            /// \param m A fwdpp::ts::marginal_tree
            /// \param n Index of a fwdpp::ts::node
            /// \param left_to_right_traversal.  If `true`, traverse children left-to-right. Else, the opposite.
            {
            }

            inline table_index_t
            operator()()
            /// Advance to the next child.
			///
            /// Returns fwdpp::ts::NULL_INDEX when no more children remain.
            {
                static_assert(NULL_INDEX < 0,
                        "NULL_INDEX < 0 is false, so something needs to change here");
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
            /// \param f A function equivalent to void (*process_child)(table_index_t)
			/// 
            /// Returns false to signify end of iteration.
            {
                auto c = this->operator()();
                bool rv = (c != NULL_INDEX);
                if (rv)
                    {
                        f(c);
                    }
                return rv;
            }
        };

        template <typename F>
        inline void
        process_children(const marginal_tree& m, table_index_t n,
                         bool left_to_right_traversal, const F& f)
		/// \brief Apply a function to children of node \a n
		/// \param m A fwdpp::ts::marginal_tree
		/// \param n Index a fwdpp::ts::node
		/// \param left_to_right_traversal.  If `true`, traverse children left-to-right. Else, the opposite.
		/// \param f A function equivalent to void (*process_child)(table_index_t)
        {
            child_iterator ci(m, n, left_to_right_traversal);
            while (ci(f))
                {
                }
        }

        inline std::vector<table_index_t>
        get_children(const marginal_tree& m, table_index_t n,
                     bool left_to_right_traversal)
		/// \brief Get all children of a node
		/// \param m A fwdpp::ts::marginal_tree
		/// \param n Index a fwdpp::ts::node
		/// \param left_to_right_traversal.  If `true`, traverse children left-to-right. Else, the opposite.
		/// \return std::vector<fwdpp::ts::table_index_t>
        {
            std::vector<table_index_t> rv;
            process_children(m, n, left_to_right_traversal,
                             [&rv](const table_index_t c) { rv.push_back(c); });
            return rv;
        }

        inline int
        num_children(const marginal_tree& m, table_index_t n)
		/// \brief Return number of children of a node
		/// \param m A fwdpp::ts::marginal_tree
		/// \param n Index a fwdpp::ts::node
		/// \return int
        {
            int nc = 0;
            process_children(m, n, true,
                             [&nc](const table_index_t /*c*/) { ++nc; });
            return nc;
        }
    } // namespace ts
} // namespace fwdpp

#endif
