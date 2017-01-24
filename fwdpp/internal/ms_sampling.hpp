#ifndef __FWDPP_INTERNAL_MS_SAMPLING_HPP__
#define __FWDPP_INTERNAL_MS_SAMPLING_HPP__

#include <cassert>

namespace KTfwd
{
    namespace fwdpp_internal
    {
        inline void
        remove_no_derived(sample_t *block)
        {
            block->erase(std::remove_if(block->begin(), block->end(),
                                        [](std::pair<double, std::string> &p) {
                                            return unsigned(std::count(
                                                       p.second.begin(),
                                                       p.second.end(), '0'))
                                                   == p.second.size();
                                        }),
                         block->end());
        }
        // Used when nsam is odd.  We just clip off the last individual
        inline void
        trim_last(sample_t *block)
        {
            std::for_each(block->begin(), block->end(),
                          [](std::pair<double, std::string> &p) {
                              if (!p.second.empty())
                                  {
                                      // remove last character
                                      p.second.erase(p.second.end() - 1);
                                  }
                          });
            remove_no_derived(block);
        }

        template <typename gamete_mcont_t, typename mcont_t,
                  typename pos_finder>
        void
        update_sample_block(sample_t &block, const gamete_mcont_t &variants,
                            const mcont_t &mutations, const size_t &i,
                            const size_t &n, const pos_finder &pf,
                            const size_t &offset = 0, const size_t &scalar = 2)
        {
            for (const auto &m : variants)
                {
                    auto mutpos = mutations[m].pos;
                    auto itr = std::find_if(
                        block.begin(), block.end(),
                        std::bind(pf, std::placeholders::_1, mutpos));
                    if (itr == block.end())
                        {
                            block.push_back(
                                std::make_pair(mutpos, std::string(n, '0')));
                            block[block.size() - 1].second[scalar * i + offset]
                                = '1';
                        }
                    else
                        {
                            itr->second[scalar * i + offset] = '1';
                        }
                }
        }

        inline void
        remove_fixed_variants_from_sample(
            std::vector<sample_site_t> &sample,
            const std::vector<sample_site_t>::size_type nsam)
        {
            sample.erase(std::remove_if(sample.begin(), sample.end(),
                                        [nsam](const sample_site_t &site) {
                                            assert(site.second.size() == nsam);
                                            return unsigned(std::count(
                                                       site.second.begin(),
                                                       site.second.end(), '1'))
                                                   == nsam;
                                        }),
                         sample.end());
        }

        template <typename mcont_t, typename gcont_t, typename dipvector_t,
                  typename integer_type = std::size_t>
        sep_sample_t
        ms_sample_separate_single_deme(
            const mcont_t &mutations, const gcont_t &gametes,
            const dipvector_t &diploids,
            const std::vector<integer_type> &diplist, const unsigned &n,
            const bool &remove_fixed)
        {
            sep_sample_t rv;
            sample_t::iterator itr;

            std::function<bool(const sample_site_t &, const double &)>
                sitefinder = [](const sample_site_t &site, const double &d) {
                    return std::fabs(site.first - d)
                           <= std::numeric_limits<double>::epsilon();
                };

            for (size_t i = 0; i < diplist.size(); ++i)
                {
                    typename dipvector_t::difference_type ind = diplist[i];
                    assert(ind >= 0);
                    assert(unsigned(ind) < diploids.size());
                    fwdpp_internal::update_sample_block(
                        rv.first, gametes[diploids[ind].first].mutations,
                        mutations, i, 2 * diplist.size(), sitefinder);
                    fwdpp_internal::update_sample_block(
                        rv.first, gametes[diploids[ind].second].mutations,
                        mutations, i, 2 * diplist.size(), sitefinder, 1);
                    fwdpp_internal::update_sample_block(
                        rv.second, gametes[diploids[ind].first].smutations,
                        mutations, i, 2 * diplist.size(), sitefinder);
                    fwdpp_internal::update_sample_block(
                        rv.second, gametes[diploids[ind].second].smutations,
                        mutations, i, 2 * diplist.size(), sitefinder, 1);
                }
            if (remove_fixed && !rv.first.empty())
                {
                    remove_fixed_variants_from_sample(rv.first,
                                                      2 * diplist.size());
                }
            if (!rv.first.empty())
                {
                    std::sort(rv.first.begin(), rv.first.end(),
                              [](const sample_site_t &lhs,
                                 const sample_site_t &rhs) {
                                  return lhs.first < rhs.first;
                              });
                }
            if (remove_fixed && !rv.second.empty())
                {
                    remove_fixed_variants_from_sample(rv.second,
                                                      2 * diplist.size());
                }
            if (!rv.second.empty())
                {
                    std::sort(rv.second.begin(), rv.second.end(),
                              [](const sample_site_t &lhs,
                                 const sample_site_t &rhs) {
                                  return lhs.first < rhs.first;
                              });
                }
            // Deal w/odd sample sizes
            if (n % 2 != 0.)
                {
                    trim_last(&rv.first);
                    trim_last(&rv.second);
                }
            assert(std::is_sorted(
                rv.first.begin(), rv.first.end(),
                [](const sample_site_t &a, const sample_site_t &b) noexcept {
                    return a.first < b.first;
                }));
            assert(std::is_sorted(
                rv.second.begin(), rv.second.end(),
                [](const sample_site_t &a, const sample_site_t &b) noexcept {
                    return a.first < b.first;
                }));
            return rv;
        }

        template <typename mcont_t, typename gcont_t, typename dipvector_t,
                  typename integer_type = std::size_t>
        std::vector<sep_sample_t>
        ms_sample_separate_mlocus(const mcont_t &mutations,
                                  const gcont_t &gametes,
                                  const dipvector_t &diploids,
                                  const std::vector<integer_type> &diplist,
                                  const unsigned &n, const bool &remove_fixed)
        {
            assert(!diploids.empty());
            using rvtype = std::vector<sep_sample_t>;
            // using genotype = typename dipvector_t::value_type;

            rvtype rv(diploids[0].size());

            std::function<bool(const sample_site_t &, const double &)>
                sitefinder = [](const sample_site_t &site, const double &d) {
                    return std::fabs(site.first - d)
                           <= std::numeric_limits<double>::epsilon();
                };

            // Go over each indidivual's mutations and update the return value
            // typename dipvector_t::const_iterator dbegin = diploids->begin();
            for (unsigned ind = 0; ind < diplist.size(); ++ind)
                {
                    unsigned rv_count = 0;
                    for (const auto &locus : diploids[diplist[ind]])
                        {
                            // finally, we can go over mutations
                            fwdpp_internal::update_sample_block(
                                rv[rv_count].first,
                                gametes[locus.first].mutations, mutations, ind,
                                2 * diplist.size(), sitefinder);
                            fwdpp_internal::update_sample_block(
                                rv[rv_count].second,
                                gametes[locus.first].smutations, mutations,
                                ind, 2 * diplist.size(), sitefinder);
                            fwdpp_internal::update_sample_block(
                                rv[rv_count].first,
                                gametes[locus.second].mutations, mutations,
                                ind, 2 * diplist.size(), sitefinder, 1);
                            fwdpp_internal::update_sample_block(
                                rv[rv_count].second,
                                gametes[locus.second].smutations, mutations,
                                ind, 2 * diplist.size(), sitefinder, 1);
                            ++rv_count;
                        }
                }

            if (remove_fixed)
                {
                    for (unsigned i = 0; i < rv.size(); ++i)
                        {
                            remove_fixed_variants_from_sample(
                                rv[i].first, 2 * diplist.size());
                            remove_fixed_variants_from_sample(
                                rv[i].second, 2 * diplist.size());
                        }
                }
            // sort on position
            for (unsigned i = 0; i < rv.size(); ++i)
                {
                    std::sort(rv[i].first.begin(), rv[i].first.end(),
                              [](const sample_site_t &lhs,
                                 const sample_site_t &rhs) {
                                  return lhs.first < rhs.first;
                              });
                    std::sort(rv[i].second.begin(), rv[i].second.end(),
                              [](const sample_site_t &lhs,
                                 const sample_site_t &rhs) {
                                  return lhs.first < rhs.first;
                              });
                    // Deal w/odd sample sizes
                    if (n % 2 != 0.)
                        {
                            trim_last(&rv[i].first);
                            trim_last(&rv[i].second);
                        }
                }
#ifndef NDEBUG
            for (const auto &rvi : rv)
                {
                    assert(std::is_sorted(rvi.first.begin(), rvi.first.end(),
                                          [](const sample_site_t &a,
                                             const sample_site_t &b) noexcept {
                                              return a.first < b.first;
                                          }));
                    assert(std::is_sorted(rvi.second.begin(), rvi.second.end(),
                                          [](const sample_site_t &a,
                                             const sample_site_t &b) noexcept {
                                              return a.first < b.first;
                                          }));
                }
#endif
            return rv;
        }
    }
}
#endif
