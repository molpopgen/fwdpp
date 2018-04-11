#ifndef __FWDPP_SUGAR_SAMPLING_HPP__
#define __FWDPP_SUGAR_SAMPLING_HPP__
#include <stdexcept>
#include <type_traits>
#include <algorithm>
#include <stdexcept>
#include <fwdpp/type_traits.hpp>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/sampling_functions.hpp>
#include <fwdpp/sugar/slocuspop.hpp>
#include <fwdpp/sugar/mlocuspop.hpp>
#include <fwdpp/sugar/sampling/sampling_details.hpp>

namespace fwdpp
{
    template <typename poptype>
    typename std::enable_if<std::is_same<typename poptype::popmodel_t,
                                         sugar::SINGLELOC_TAG>::value,
                            sample_t>::type
    sample(const gsl_rng *r, const poptype &p, const unsigned nsam,
           const bool removeFixed)
    /*!
      Take a random sample of nsam chromosomes from a population

      \param r A random-number generator
      \param p A population
      \param nsam The sample size
      \param removeFixed Whether or not to remove variants present in all nsam
      chromosomes

      \return A vector of both neutral and non-neutral variants
    */
    {
        static_assert(std::is_same<typename poptype::popmodel_t,
                                   sugar::SINGLELOC_TAG>::value,
                      "poptype must be SINGLELOC_TAG");
        auto rv = ms_sample(r, p.mutations, p.gametes, p.diploids, nsam,
                            removeFixed);
        finish_sample(rv, p.fixations, nsam, removeFixed,
                      sugar::treat_neutral::ALL);
        return rv;
    }

    template <typename poptype>
    typename std::enable_if<std::is_same<typename poptype::popmodel_t,
                                         sugar::MULTILOC_TAG>::value,
                            std::vector<sample_t>>::type
    sample(const gsl_rng *r, const poptype &p, const unsigned nsam,
           const bool removeFixed)
    /*!
      Take a random sample of nsam chromosomes from a population

      \param r A random-number generator
      \param p A population
      \param nsam The sample size
      \param removeFixed Whether or not to remove variants present in all nsam
      chromosomes

      \return A vector of both neutral and non-neutral variants
    */
    {
        static_assert(std::is_same<typename poptype::popmodel_t,
                                   sugar::MULTILOC_TAG>::value,
                      "poptype must be MULTILOC_TAG");
        if (!removeFixed && p.locus_boundaries.empty())
            {
                throw std::invalid_argument(
                    "locus boundaries required when adding fixations");
            }
        auto rv = ms_sample(r, p.mutations, p.gametes, p.diploids, nsam,
                            removeFixed);
        if (!removeFixed && p.locus_boundaries.size() != rv.size())
            {
                throw std::invalid_argument(
                    "incorrect number of elements in locus_boundaries");
            }
        finish_sample(rv, p.fixations, nsam, removeFixed,
                      sugar::treat_neutral::ALL, p.locus_boundaries);
        return rv;
    }

    template <typename poptype>
    typename std::enable_if<std::is_same<typename poptype::popmodel_t,
                                         sugar::MULTILOC_TAG>::value,
                            std::vector<sep_sample_t>>::type
    sample_separate(
        const gsl_rng *r, const poptype &p, const unsigned nsam,
        const bool removeFixed)
    /*!
      Take a random sample of nsam chromosomes from a population

      \param r A random-number generator
      \param p A population
      \param nsam The sample size
      \param removeFixed Whether or not to remove variants present in all nsam
      chromosomes

      \return A pair of vectors.  The first element contains neutral variants.
      The second contains non-neutral variants.
    */
    {
        static_assert(std::is_same<typename poptype::popmodel_t,
                                   sugar::MULTILOC_TAG>::value,
                      "poptype must be MULTILOC_TAG");
        if (!removeFixed && p.locus_boundaries.empty())
            {
                throw std::invalid_argument(
                    "locus boundaries required when adding fixations");
            }
        auto rv = ms_sample_separate(r, p.mutations, p.gametes, p.diploids,
                                     nsam, removeFixed);
        if (!removeFixed && p.locus_boundaries.size() != rv.size())
            {
                throw std::invalid_argument(
                    "incorrect number of elements in locus_boundaries");
            }
        finish_sample(rv, p.fixations, nsam, removeFixed,
                      sugar::treat_neutral::ALL, p.locus_boundaries);
        return rv;
    }

    template <typename poptype>
    typename std::enable_if<std::is_same<typename poptype::popmodel_t,
                                         sugar::SINGLELOC_TAG>::value,
                            sep_sample_t>::type
    sample_separate(const gsl_rng *r, const poptype &p, const unsigned nsam,
                    const bool removeFixed)
    /*!
      Take a random sample of nsam chromosomes from a population

      \param r A random-number generator
      \param p A population
      \param nsam The sample size
      \param removeFixed Whether or not to remove variants present in all nsam
      chromosomes

      \return A pair of vectors.  The first element contains neutral variants.
      The second contains non-neutral variants.
    */
    {
        static_assert(std::is_same<typename poptype::popmodel_t,
                                   sugar::SINGLELOC_TAG>::value,
                      "poptype must be SINGLELOC_TAG");
        auto rv = ms_sample_separate(r, p.mutations, p.gametes, p.diploids,
                                     nsam, removeFixed);
        finish_sample(rv, p.fixations, nsam, removeFixed,
                      sugar::treat_neutral::ALL);
        return rv;
    }

    template <typename poptype, typename integer_type = std::size_t>
    typename std::enable_if<std::is_same<typename poptype::popmodel_t,
                                         sugar::SINGLELOC_TAG>::value,
                            sample_t>::type
    sample(const poptype &p, const std::vector<integer_type> &individuals,
           const bool removeFixed)
    /*!
      Take a non-random sample of diploids from a population

      \param p A population
      \param individuals The indexes of the diploids to sample
      \param removeFixed Whether or not to remove variants present in all nsam
      chromosomes

      \return A vector of both neutral and non-neutral variants
    */
    {
        static_assert(std::is_same<typename poptype::popmodel_t,
                                   sugar::SINGLELOC_TAG>::value,
                      "poptype must be SINGLELOC_TAG");
        if (individuals.empty())
            return sample_t();
        if (std::find_if(
                individuals.begin(), individuals.end(),
                [&p](const unsigned &u) { return u >= p.diploids.size(); })
            != individuals.end())
            {
                throw std::out_of_range(
                    "fwdpp::sample_separate: individual index out of range");
            }

        auto rv = sample_details(p, individuals, removeFixed);
        return rv;
    }

    template <typename poptype, typename integer_type = std::size_t>
    typename std::enable_if<std::is_same<typename poptype::popmodel_t,
                                         sugar::MULTILOC_TAG>::value,
                            std::vector<sample_t>>::type
    sample(const poptype &p, const std::vector<integer_type> &individuals,
           const bool removeFixed)
    /*!
      Take a non-random sample of diploids from a population

      \param p A population
      \param individuals The indexes of the diploids to sample
      \param removeFixed Whether or not to remove variants present in all nsam
      chromosomes

      \return A vector of both neutral and non-neutral variants
    */
    {
        static_assert(std::is_same<typename poptype::popmodel_t,
                                   sugar::MULTILOC_TAG>::value,
                      "poptype must be MULTILOC_TAG");
        if (!removeFixed && p.locus_boundaries.empty())
            {
                throw std::invalid_argument(
                    "locus boundaries required when adding fixations");
            }
        if (individuals.empty())
            return std::vector<sample_t>();
        if (std::find_if(
                individuals.begin(), individuals.end(),
                [&p](const unsigned &u) { return u >= p.diploids.size(); })
            != individuals.end())
            {
                throw std::out_of_range(
                    "fwdpp::sample_separate: individual index out of range");
            }

        return sample_details(p, individuals, removeFixed, p.locus_boundaries);
    }

    template <typename poptype, typename integer_type = std::size_t>
    typename std::enable_if<std::is_same<typename poptype::popmodel_t,
                                         sugar::SINGLELOC_TAG>::value,
                            sep_sample_t>::type
    sample_separate(const poptype &p,
                    const std::vector<integer_type> &individuals,
                    const bool removeFixed)
    /*!
      Take a non-random sample of nsam chromosomes from a population

      \param p A population
      \param individuals The indexes of the diploids to 'view' from the
      population
      \param removeFixed Whether or not to remove variants present in all nsam
      chromosomes

      \return A pair of vectors.  The first element contains neutral variants.
      The second contains non-neutral variants.
    */
    {
        static_assert(std::is_same<typename poptype::popmodel_t,
                                   sugar::SINGLELOC_TAG>::value,
                      "poptype must be SINGLELOC_TAG");
        if (individuals.empty())
            return sep_sample_t();
        if (std::find_if(
                individuals.begin(), individuals.end(),
                [&p](const unsigned &u) { return u >= p.diploids.size(); })
            != individuals.end())
            {
                throw std::out_of_range(
                    "fwdpp::sample_separate: individual index out of range");
            }
        auto rv = fwdpp_internal::ms_sample_separate_single_locus_pop(
            p.mutations, p.gametes, p.diploids, individuals,
            2 * individuals.size(), removeFixed);
        finish_sample(rv, p.fixations, 2 * individuals.size(), removeFixed,
                      sugar::treat_neutral::ALL);
        return rv;
    }

    template <typename poptype, typename integer_type = std::size_t>
    typename std::enable_if<std::is_same<typename poptype::popmodel_t,
                                         sugar::MULTILOC_TAG>::value,
                            std::vector<sep_sample_t>>::type
    sample_separate(
        const poptype &p, const std::vector<integer_type> &individuals,
        const bool removeFixed)
    /*!
      Take a non-random sample of nsam chromosomes from a population

      \param p A population
      \param individuals The indexes of the diploids to 'view' from the
      population
      \param removeFixed Whether or not to remove variants present in all nsam
      chromosomes

      \return A pair of vectors.  The first element contains neutral variants.
      The second contains non-neutral variants.
    */
    {
        static_assert(std::is_same<typename poptype::popmodel_t,
                                   sugar::MULTILOC_TAG>::value,
                      "poptype must be MULTILOC_TAG");
        if (!removeFixed && p.locus_boundaries.empty())
            {
                throw std::invalid_argument(
                    "locus boundaries required when adding fixations");
            }
        if (individuals.empty())
            return std::vector<sep_sample_t>();
        if (std::find_if(
                individuals.begin(), individuals.end(),
                [&p](const unsigned &u) { return u >= p.diploids.size(); })
            != individuals.end())
            {
                throw std::out_of_range(
                    "fwdpp::sample_separate: individual index out of range");
            }
        auto rv = fwdpp_internal::ms_sample_separate_mlocus(
            p.mutations, p.gametes, p.diploids, individuals,
            2 * individuals.size(), removeFixed);
        if (!removeFixed && p.locus_boundaries.size() != rv.size())
            {
                throw std::invalid_argument(
                    "incorrect number of elements in locus_boundaries");
            }
        finish_sample(
            rv, p.fixations, 2 * static_cast<unsigned>(individuals.size()),
            removeFixed, sugar::treat_neutral::ALL, p.locus_boundaries);
        return rv;
    }
}

#endif
