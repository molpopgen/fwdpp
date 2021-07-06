#ifndef FWDPP_TS_MARGINAL_TREE_FUNCTIONS_ROOTS_ITERATOR_HPP
#define FWDPP_TS_MARGINAL_TREE_FUNCTIONS_ROOTS_ITERATOR_HPP

#include <type_traits>
#include <stdexcept>
#include "../marginal_tree.hpp"
#include "../types/generate_null_id.hpp"

namespace fwdpp
{
    namespace ts
    {
        template <typename SignedInteger> class root_iterator
        /// Allows iteration over tree roots
        /// \note Tied to lifetime of a marginal_tree.
        /// \headerfile fwdpp/ts/marginal_tree_functions/roots.hpp
        {
          private:
            SignedInteger current_root;
            typename std::vector<SignedInteger>::const_iterator rsib_beg, rsib_end;

          public:
            explicit root_iterator(const marginal_tree<SignedInteger>& m)
                : current_root(m.left_root), rsib_beg(begin(m.right_sib)),
                  rsib_end(end(m.right_sib))
            /// \param m A fwdpp::ts::marginal_tree
            {
                if (current_root == types::generate_null_id<SignedInteger>())
                    {
                        throw std::invalid_argument("root list is empty");
                    }
                if (rsib_beg == rsib_end)
                    {
                        throw std::invalid_argument("empty list of right sibs");
                    }
            }

            inline SignedInteger
            operator()()
            /// \return The next root notde
            /// NULL_INDEX signals end of iteration
            {
                static_assert(
                    types::generate_null_id<SignedInteger>() < 0,
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
                bool rv = (croot != types::generate_null_id<SignedInteger>());
                if (rv)
                    {
                        f(croot);
                    }
                return rv;
            }
        };

        template <typename SignedInteger, typename F>
        inline void
        process_roots(const marginal_tree<SignedInteger>& m, const F& f)
        /// \brief Apply a function to all roots
        /// \param m A Fwdpp::ts::marginal_tree
        /// \param f A function equivalent to void(*process_root)(SignedInteger)
        {
            root_iterator<SignedInteger> ri(m);
            while (ri(f))
                {
                }
        }

        template <typename SignedInteger>
        inline std::vector<SignedInteger>
        get_roots(const marginal_tree<SignedInteger>& m)
        /// \brief Get a vector of all roots
        /// \param m A Fwdpp::ts::marginal_tree
        {
            std::vector<SignedInteger> rv;
            process_roots(m, [&rv](const SignedInteger r) { rv.push_back(r); });
            return rv;
        }

        template <typename SignedInteger>
        inline int
        num_roots(const marginal_tree<SignedInteger>& m)
        /// \brief Return number of roots
        /// \param m A Fwdpp::ts::marginal_tree
        {
            int nroots = 0;
            process_roots(m, [&nroots](const SignedInteger /*r*/) { ++nroots; });
            return nroots;
        }

    } // namespace ts
} // namespace fwdpp
#endif
