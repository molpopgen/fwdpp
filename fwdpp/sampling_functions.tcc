//  -*- C++ -*-
#ifndef __FWDPP_SAMPLING_FUNCTIONS_TCC__
#define __FWDPP_SAMPLING_FUNCTIONS_TCC__

#include <fwdpp/fwd_functional.hpp>
#include <fwdpp/internal/ms_sampling.hpp>
#include <limits>
#include <type_traits>
#include <algorithm>

namespace KTfwd
{
    template <typename gamete_type, typename allocator_t,
              template <typename, typename> class container_type>
    std::vector<unsigned>
    sample(const gsl_rng *r,
           const container_type<gamete_type, allocator_t> &gametes,
           const unsigned &n, const unsigned &N)
    {
        std::vector<double> freqs;
        std::vector<unsigned> counts(gametes.size(), 0);
        std::for_each(gametes.begin(), gametes.end(),
                      [&freqs, &N](const gamete_type &__g) {
                          freqs.emplace_back(
                              std::move(double(__g.n) / double(N)));
                      });
        gsl_ran_multinomial(r, gametes.size(), n, &freqs[0], &counts[0]);
        return counts;
    }

    // SAMPLERS FOR INDIVIDUAL-BASED SIMULATIONS
    template <typename mcont_t, typename gcont_t, typename allocator,
              typename diploid_geno_t,
              template <typename, typename> class vector_type>
    typename std::enable_if<traits::is_diploid<diploid_geno_t>::value,
                            sample_t>::type
    ms_sample(const gsl_rng *r, const mcont_t &mutations,
              const gcont_t &gametes,
              const vector_type<diploid_geno_t, allocator> &diploids,
              const unsigned &n, const bool &remove_fixed)
    {
        auto separate = ms_sample_separate(r, mutations, gametes, diploids, n,
                                           remove_fixed);
        std::move(separate.second.begin(), separate.second.end(),
                  std::back_inserter(separate.first));
        std::sort(separate.first.begin(), separate.first.end(),
                  [](const sample_site_t &lhs, const sample_site_t &rhs) {
                      return lhs.first < rhs.first;
                  });
        return separate.first;
    }

    template <typename mcont_t, typename gcont_t, typename allocator,
              typename diploid_geno_t,
              template <typename, typename> class vector_type>
    typename std::enable_if<traits::is_diploid<diploid_geno_t>::value,
                            sep_sample_t>::type
    ms_sample_separate(const gsl_rng *r, const mcont_t &mutations,
                       const gcont_t &gametes,
                       const vector_type<diploid_geno_t, allocator> &diploids,
                       const unsigned &n, const bool &remove_fixed)
    {
        std::vector<unsigned> diplist;
        unsigned isodd = (n % 2 != 0.) ? 1u : 0u;
        for (unsigned i = 0; i < n / 2 + isodd; ++i)
            {
                diplist.push_back(std::vector<unsigned>::value_type(
                    gsl_ran_flat(r, 0., double(diploids.size()))));
            }
        return fwdpp_internal::ms_sample_separate_single_deme(
            mutations, gametes, diploids, diplist, n, remove_fixed);
    }

    // Individual-based sims, multilocus algorithm
    template <typename mcont_t, typename gcont_t, typename dcont_t>
    typename std::enable_if<traits::is_diploid<typename dcont_t::value_type::
                                                   value_type>::value,
                            std::vector<sep_sample_t>>::type
    ms_sample_separate(const gsl_rng *r, const mcont_t &mutations,
                       const gcont_t &gametes, const dcont_t &diploids,
                       const unsigned &n, const bool &remove_fixed)
    {
        std::vector<unsigned> diplist;
        unsigned isodd = (n % 2 != 0.) ? 1u : 0u;
        for (unsigned i = 0; i < n / 2 + isodd; ++i)
            {
                diplist.push_back(
                    unsigned(gsl_ran_flat(r, 0., double(diploids.size()))));
            }
        return fwdpp_internal::ms_sample_separate_mlocus(
            mutations, gametes, diploids, diplist, n, remove_fixed);
    }

    template <typename mcont_t, typename gcont_t, typename dcont_t>
    typename std::enable_if<traits::is_diploid<typename dcont_t::value_type::
                                                   value_type>::value,
                            std::vector<sample_t>>::type
    ms_sample(const gsl_rng *r, const mcont_t &mutations,
              const gcont_t &gametes, const dcont_t &diploids,
              const unsigned &n, const bool &remove_fixed)
    {
        auto separate = ms_sample_separate(r, mutations, gametes, diploids, n,
                                           remove_fixed);
        std::vector<sample_t> rv;
        for (unsigned i = 0; i < separate.size(); ++i)
            {
                std::move(separate[i].second.begin(), separate[i].second.end(),
                          std::back_inserter(separate[i].first));
                std::sort(
                    separate[i].first.begin(), separate[i].first.end(),
                    [](const sample_site_t &lhs, const sample_site_t &rhs) {
                        return lhs.first < rhs.first;
                    });
                rv.emplace_back(std::move(separate[i].first));
            }
        return rv;
    }
}

#endif
