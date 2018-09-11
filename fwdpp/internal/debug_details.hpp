#ifndef FWDPP_INTERNAL_DEBUG_DETAILS_HPP
#define FWDPP_INTERNAL_DEBUG_DETAILS_HPP

#include <algorithm>
#include <stdexcept>

namespace fwdpp
{
    namespace debug
    {
        namespace detail
        {
#ifndef NDEBUG
            template <typename mcont_t, typename iterator>
            void
            validate_mutation_key_ranges(const mcont_t &mutations,
                                         const iterator beg,
                                         const iterator end)
            {
                if (std::any_of(beg, end, [&mutations](const std::size_t m) {
                        return m >= mutations.size();
                    }))
                    {
                        throw std::runtime_error("FWDPP DEBUG: mutation key "
                                                 "out of range");
                    }
            }

            template <typename gamete_t, typename mcont_t>
            void
            gamete_is_sorted(const gamete_t &g, const mcont_t &m)
            /*!
      \brief Check that neutral mutation keys are sorted according to mutation
      position
    */
            {
                const auto comp = [&m](const size_t i, const size_t j) {
                    return m[i].pos <= m[j].pos;
                };

                if (!std::is_sorted(g.mutations.begin(), g.mutations.end(),
                                    comp))
                    {
                        throw std::runtime_error(
                            "neutral mutation keys not sorted");
                    }
                if (!std::is_sorted(g.smutations.begin(), g.smutations.end(),
                                    comp))
                    {
                        throw std::runtime_error(
                            "selected mutation keys not sorted");
                    }
            }
#else
            template <typename mcont_t, typename iterator>
            void
            validate_mutation_key_ranges(const mcont_t & /*mutations*/,
                                         const iterator /*beg*/,
                                         const iterator /*end*/)
            {
            }

            template <typename gamete_t, typename mcont_t>
            void
            gamete_is_sorted(const gamete_t &, const mcont_t &)
            /*!
      \brief Check that neutral mutation keys are sorted according to mutation
      position
    */
            {
            }
#endif
        }
    }
}

#endif
