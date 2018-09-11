#ifndef FWDPP_DEBUG_HPP
#define FWDPP_DEBUG_HPP

#include <algorithm>
#include <numeric>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/type_traits.hpp>
#include "internal/debug_details.hpp"

namespace fwdpp
{
    namespace debug
    {
        template <typename gcont_t>
        void validate_sum_gamete_counts(const gcont_t &gametes,
                                        const uint_t expected_sum)
        {
            detail::validate_sum_gamete_counts(gametes, expected_sum);
        }

        template <typename mcont_t, typename iterator>
        void
        validate_mutation_key_ranges(const mcont_t &mutations,
                                     const iterator beg, const iterator end)
        /*! Throw an exception if any mutation keys are >= mutations.size()
         */
        {
            detail::validate_mutation_key_ranges(mutations, beg, end);
        }

        template <typename gamete_t, typename mcont_t>
        void
        gamete_is_sorted(const gamete_t &g, const mcont_t &mutations)
        /*!
          \brief Check that neutral mutation keys are sorted according to mutation
          position
        */
        {
            detail::gamete_is_sorted(g, mutations);
        }

        template <typename gamete_t, typename mcont_t>
        void
        gamete_data_valid(const gamete_t &g, const mcont_t &mutations,
                          const std::vector<uint_t> &mutcounts)
        {
            detail::gamete_data_valid(g, mutations, mutcounts);
        }

        template <typename poptype>
        void
        validate_pop_data(const poptype &pop)
        {
            detail::validate_pop_data(pop);
        }
    } // namespace debug
} // namespace fwdpp

#endif
