#ifndef FWDPP_DEBUG_HPP
#define FWDPP_DEBUG_HPP

#include <algorithm>
#include <numeric>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/type_traits.hpp>
#include "internal/debug_details.hpp"

/*! \namespace fwdpp::debug Debugging routines
 * In general, the library does not throw exceptions
 * during the execution of simulation-related algorithms.
 * Rather, exception handling is limited to sanity-checking
 * parameters inputs, etc., for object construction and other routines.
 *
 * However, it is sometimes useful to perform expensive checks on the
 * data integrity of a simulated population.  The routines in
 * fwdpp::debug allow you to do that without affecting the 
 * performance of code compiled in "release" mode.
 *
 * Following the tradition of the C language, compiling with
 * -DNDEBUG signals that a release-mode build is desired, 
 * meaning that expensive runtime checks are disables.
 * 
 * If that symbol is not defined, then various runtime checks
 * are enabled, and any failed checks with throw std::runtime_error.
 * The message field of the error will begin with "FWDPP DEBUG:".
 *
 * Functions in namespace fwdpp::debug will be optimized out when
 * code is compiled with -DNDEBUG, as they evaluate to empty functions.
 */

namespace fwdpp
{
    namespace debug
    {
        template <typename gcont_t>
        void
        validate_sum_gamete_counts(const gcont_t &gametes,
                                   const uint_t expected_sum)
        /// Throw exception if sum of all gamete counts != \a expected_sum
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

        template <typename gamete_t>
        void
        gamete_is_extant(const gamete_t &gamete)
        /// Throw exception of gamete.n == 0
        {
            detail::gamete_is_extant(gamete);
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
        /// Throw exception if gamete data are unsorted or if mutation
        /// key storage is invalid.
        /// \note Only call this if the gamete is extant!!!!
        {
            detail::gamete_data_valid(g, mutations, mutcounts);
        }

        template <typename poptype>
        void
        validate_pop_data(const poptype &pop)
        /// Throw exception if any in a series of checks on
        /// the internal state of \a pop fail.
        /// \note This is an expensive function!
        {
            detail::validate_pop_data(pop);
        }

        template <typename mutation_type>
        void
        check_mutation_neutrality(const mutation_type &mutation,
                                  const bool expected_neutrality)
        /// Throw exception if mutation.neutral != \a expected_neutrality
        {
            detail::check_mutation_neutrality(mutation, expected_neutrality);
        }

        template <typename poptype>
        /// Throw exception if gametes contain keys to extinct gametes
        /// \version 0.7.4 Added to fwdpp
        void
        all_gametes_extant(const poptype &pop)
        {
            detail::all_gametes_extant(pop, typename poptype::popmodel_t());
        }
    } // namespace debug
} // namespace fwdpp

#endif
