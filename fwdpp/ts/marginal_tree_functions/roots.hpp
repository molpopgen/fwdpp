#ifndef FWDPP_TS_MARGINAL_TREE_FUNCTIONS_ROOTS_ITERATOR_HPP
#define FWDPP_TS_MARGINAL_TREE_FUNCTIONS_ROOTS_ITERATOR_HPP

#include <type_traits>
#include <stdexcept>
#include "../marginal_tree.hpp"

namespace fwdpp
{
    namespace ts
    {
        class root_iterator
        /// Allows iteration over tree roots
        /// \note Tied to lifetime of a marginal_tree.
        /// \headerfile fwdpp/ts/marginal_tree_functions/roots.hpp
        {
          private:
            table_index_t current_root;
            std::vector<table_index_t>::const_iterator rsib_beg, rsib_end;

          public:
            explicit root_iterator(const marginal_tree& m)
                : current_root(m.left_root), rsib_beg(begin(m.right_sib)),
                  rsib_end(end(m.right_sib))
            /// \param m A fwdpp::ts::marginal_tree
            {
                if (current_root == NULL_INDEX)
                    {
                        throw std::invalid_argument("root list is empty");
                    }
                if (rsib_beg == rsib_end)
                    {
                        throw std::invalid_argument(
                            "empty list of right sibs");
                    }
            }

            inline table_index_t
            operator()()
            /// \return The next root notde
            /// NULL_INDEX signals end of iteration
            {
                static_assert(NULL_INDEX < 0,
                        "NULL_INDEX < 0 is false, so something needs to change here");
                auto croot = current_root;
                if (croot < 0)
                    {
                        return croot;
                    }
                if (rsib_beg + croot >= rsib_end)
                    {
                        throw std::runtime_error("root iteration error");
                    }
                current_root = *(rsib_beg + croot);
                return croot;
            }

            template <typename F>
            inline bool
            operator()(const F& f)
            /// Apply a function to all roots
            /// \return true if there are more roots to iterate over, false otherwise.
            {
                auto croot = this->operator()();
                bool rv = (croot != NULL_INDEX);
                if (rv)
                    {
                        f(croot);
                    }
                return rv;
            }
        };

        template <typename F>
        inline void
        process_roots(const marginal_tree& m, const F& f)
        /// \brief Apply a function to all roots
        /// \param m A Fwdpp::ts::marginal_tree
        /// \param f A function equivalent to void(*process_root)(table_index_t)
        {
            root_iterator ri(m);
            while (ri(f))
                {
                }
        }

        inline std::vector<table_index_t>
        get_roots(const marginal_tree& m)
        /// \brief Get a vector of all roots
        /// \param m A Fwdpp::ts::marginal_tree
        {
            std::vector<table_index_t> rv;
            process_roots(m, [&rv](const table_index_t r) { rv.push_back(r); });
            return rv;
        }

        inline int
        num_roots(const marginal_tree& m)
        /// \brief Return number of roots
        /// \param m A Fwdpp::ts::marginal_tree
        {
            int nroots = 0;
            process_roots(m, [&nroots](const table_index_t /*r*/) { ++nroots; });
            return nroots;
        }

    } // namespace ts
} // namespace fwdpp
#endif
