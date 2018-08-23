#ifndef __SAMPLING_FUNCTIONS_HPP__
#define __SAMPLING_FUNCTIONS_HPP__

#include <vector>
#include <map>
#include <functional>
#include <algorithm>
#include <utility>
#include <string>
#include <cassert>
#include <type_traits>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fwdpp/type_traits.hpp>
#include <fwdpp/data_matrix.hpp>

/*! @defgroup samplingPops Functions related to taking samples from simulated
  populations
*/

namespace fwdpp
{
    template <typename poptype>
    data_matrix
    sample_individuals(const poptype &pop,
                       const std::vector<std::size_t> &individuals,
                       const bool include_neutral, const bool include_selected,
                       const bool remove_fixed)
    {
        auto keys = mutation_keys(pop, individuals, include_neutral,
                                  include_neutral);
        if (remove_fixed)
            {
                keys.first.erase(
                    std::remove_if(
                        keys.first.begin(), keys.first.end(),
                        [&individuals](
                            const std::pair<std::size_t, uint_t> &p) {
                            return p.second == 2 * individuals.size();
                        }),
                    keys.first.end());
                keys.second.erase(
                    std::remove_if(
                        keys.second.begin(), keys.second.end(),
                        [&individuals](
                            const std::pair<std::size_t, uint_t> &p) {
                            return p.second == 2 * individuals.size();
                        }),
                    keys.second.end());
            }
        const auto comp = [&pop](const std::pair<std::size_t, uint_t> &a,
                                 const std::pair<std::size_t, uint_t> &b) {
            return pop.mutations[a.first].pos < pop.mutations[b.first].pos;
        };
        std::sort(keys.first.begin(), keys.first.end(), comp);
        std::sort(keys.second.begin(), keys.second.end(), comp);
        return haplotype_matrix(pop, individuals, keys.first, keys.second);
    }

    /*!
      A variable site in a sample is a pair (pos,genotypes).

      This is equivalent to libsequence's Sequence::polymorphicSite
    */
    using sample_site_t = std::pair<double, std::string>;
    /*!
      A sample is a vector of variable sites.

      This is equivalent to libsequence's Sequence::polySiteVector.

      For this type, 'neutral' and 'selected' variants are intermingled.
    */
    using sample_t = std::vector<sample_site_t>;
    /*!
      A sample where 'neutral' and 'selected' variants are separated.

      'first' contains the 'neutral' variants, and 'second' contains the
      'selected' variants.
    */
    using sep_sample_t = std::pair<sample_t, sample_t>;

    /*!
      \brief Sampling from a population in an individual-based simulation
      \return A vector of variable sites
      \ingroup samplingPopsInd
    */
    template <typename mcont_t, typename gcont_t, typename allocator,
              typename diploid_geno_t,
              template <typename, typename> class vector_type>
    typename std::enable_if<traits::is_diploid<diploid_geno_t>::value,
                            sample_t>::type
    ms_sample(const gsl_rng *r, const mcont_t &mutations,
              const gcont_t &gametes,
              const vector_type<diploid_geno_t, allocator> &diploids,
              const unsigned &n, const bool &remove_fixed = true);

    /*!
      \brief Sampling from a population in an individual-based simulation.
      Selected and neutral mutations returned separately
      \return A pair of vectors of variable sites.  The first block is neutral
      variants, the second is non-neutral variants
      \ingroup samplingPopsInd
    */
    template <typename mcont_t, typename gcont_t, typename allocator,
              typename diploid_geno_t,
              template <typename, typename> class vector_type>
    typename std::enable_if<traits::is_diploid<diploid_geno_t>::value,
                            sep_sample_t>::type
    ms_sample_separate(const gsl_rng *r, const mcont_t &mutations,
                       const gcont_t &gametes,
                       const vector_type<diploid_geno_t, allocator> &diploids,
                       const unsigned &n, const bool &remove_fixed = true);

    /*!
      \brief Sample from an individual-based, multi-locus simulation.
      \return A vector of vectors of variable sites.  There is 1 vector per
      locus.
      \note Neutral + selected mutations intermixed
      \ingroup samplingPopsInd
    */
    template <typename mcont_t, typename gcont_t, typename dcont_t>
    typename std::enable_if<traits::is_diploid<typename dcont_t::value_type::
                                                   value_type>::value,
                            std::vector<sample_t>>::type
    ms_sample(const gsl_rng *r, const mcont_t &mutations,
              const gcont_t &gametes, const dcont_t &diploids,
              const unsigned &n, const bool &remove_fixed);

    /*!
      \brief Sample from an individual-based, multi-locus simulation.
      \return A vector of pairs of vectors of variable sites.  There is 1
      vector per locus.
      \note For each locus, the first member of the pair corresponds to neutral
      sites, the second to selected.
      \ingroup samplingPopsInd
    */
    template <typename mcont_t, typename gcont_t, typename dcont_t>
    typename std::enable_if<traits::is_diploid<typename dcont_t::value_type::
                                                   value_type>::value,
                            std::vector<sep_sample_t>>::type
    ms_sample_separate(const gsl_rng *r, const mcont_t &mutations,
                       const gcont_t &gametes, const dcont_t &diploids,
                       const unsigned &n, const bool &remove_fixed);
}
#endif
#include <fwdpp/sampling_functions.tcc>
