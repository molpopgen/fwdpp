#ifndef FWDPP_INTERNAL_DEMOGRAPHY_DETAILS
#define FWDPP_INTERNAL_DEMOGRAPHY_DETAILS

#include <vector>
#include <algorithm>
#include <cassert>
#include <gsl/gsl_randist.h>
#include <fwdpp/internal/sample_diploid_helpers.hpp>

namespace KTfwd
{
    namespace fwdpp_internal
    {
        inline std::vector<std::size_t>
        sample_individuals(const gsl_rng *r, const std::size_t N,
                           const uint_t N2, const bool with_replacement)
        {
            std::vector<std::size_t> rv;
            for (std::size_t i = 0; i < N2; ++i)
                {
                    std::size_t ind
                        = std::size_t(gsl_ran_flat(r, 0., double(N)));
                    if (!with_replacement)
                        {
                            while (std::find(rv.begin(), rv.end(), ind)
                                   != rv.end())
                                {
                                    ind = std::size_t(
                                        gsl_ran_flat(r, 0., double(N)));
                                }
                        }
                    rv.push_back(ind);
                }
            // sort in descending order
            std::sort(rv.begin(), rv.end(),
                      [](size_t i, size_t j) noexcept { return i > j; });
            return rv;
        }

        template <typename mcont_t, typename mcount_t, typename gcont_t,
                  typename vdipvector_t>
        int
        split_deme_replacement(const gsl_rng *r, const mcont_t &mutations,
                               mcount_t &mcounts, gcont_t &gametes,
                               vdipvector_t &diploids, const size_t i,
                               const uint_t N_new)
        {
            auto indlist_n
                = sample_individuals(r, diploids[i].size(), N_new, true);
            auto indlist_c = sample_individuals(
                r, diploids[i].size(), diploids[i].size() - std::size_t(N_new),
                true);

            // get rid of deme i's contribution to gamete counts
            for (const auto &dip : diploids[i])
                {
                    gametes[dip.first].n--;
                    gametes[dip.second].n--;
                }

            // copy diploids from deme i
            typename vdipvector_t::value_type deme_i(diploids[i]);

            // delete diploids from deme i
            std::size_t rsize = diploids[i].size() - std::size_t(N_new);
            diploids[i].clear();
            diploids[i].reserve(rsize);

            // add a new deme
            diploids.emplace_back(typename vdipvector_t::value_type());
            diploids[i].reserve(N_new);

            // get references to current and new deme
            auto &curr_deme = diploids[i];
            auto &new_deme = diploids[diploids.size() - 1];

            // populate new deme based on copy of deme i
            for (const auto &i : indlist_n)
                {
                    new_deme.push_back(deme_i[i]);
                }

            // populate current deme from its copy
            for (const auto &i : indlist_c)
                {
                    curr_deme.push_back(deme_i[i]);
                }

            // Need to update counts
            for (const auto &dip : new_deme)
                {
                    gametes[dip.first].n++;
                    gametes[dip.second].n++;
                }
            for (const auto &dip : curr_deme)
                {
                    gametes[dip.first].n++;
                    gametes[dip.second].n++;
                }
            fwdpp_internal::process_gametes(gametes, mutations, mcounts);
            return 0;
        }

        template <typename mcont_t, typename mcount_t, typename gcont_t,
                  typename vdipvector_t>
        int
        split_deme_no_replacement(const gsl_rng *r, const mcont_t &,
                                  mcount_t &, gcont_t &,
                                  vdipvector_t &diploids, const size_t i,
                                  const uint_t N_new)
        {
            auto indlist
                = sample_individuals(r, diploids[i].size(), N_new, false);
            // add an empty new deme
            diploids.emplace_back(typename vdipvector_t::value_type());
            // get reference to new deme and current deme
            auto &new_deme = diploids[diploids.size() - 1];
            auto &curr_deme = diploids[i];
#ifndef NDEBUG
            auto curr_deme_size = curr_deme.size();
#endif
            new_deme.reserve(N_new);
            /*
              copy individuals from deme 1 to deme 2 and erase from deme 1
              std::move won't necessarily do the trick here...
            */
            for (const auto &ind : indlist)
                {
                    new_deme.emplace_back(*(curr_deme.begin() + ind));
                    curr_deme.erase(curr_deme.begin() + ind);
                }
            assert(new_deme.size() == N_new);
            assert(curr_deme.size() == curr_deme_size - N_new);
            return 0;
        }

        template <typename mcont_t, typename mcount_t, typename gcont_t,
                  typename vdipvector_t>
        int
        split_deme_details(const gsl_rng *r, const mcont_t &mutations,
                           mcount_t &mcounts, gcont_t &gametes,
                           vdipvector_t &diploids, const size_t i,
                           const uint_t N_new, const bool replacement = false)
        {
            return (replacement)
                       ? split_deme_replacement(r, mutations, mcounts, gametes,
                                                diploids, i, N_new)
                       : split_deme_no_replacement(r, mutations, mcounts,
                                                   gametes, diploids, i,
                                                   N_new);
        }
    }
}

#endif
