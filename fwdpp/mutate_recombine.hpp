#ifndef FWDPP_MUTATE_RECOMBINE_HPP__
#define FWDPP_MUTATE_RECOMBINE_HPP__

#include <vector>
#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/internal/mutation_internal.hpp>

namespace KTfwd
{
    template <typename recombination_policy, typename gcont_t,
              typename mcont_t>
    std::vector<double>
    generate_breakpoints(const std::size_t g1, const std::size_t g2,
                         const gcont_t &gametes, const mcont_t &mutations,
                         const recombination_policy rec_pol)
    {
        auto nm1
            = gametes[g1].mutations.size() + gametes[g1].smutations.size();
        auto nm2
            = gametes[g2].mutations.size() + gametes[g2].smutations.size();
        if ((std::min(nm1, nm2) == 0 && std::max(nm1, nm2) == 1)
            || gametes[g1] == gametes[g2])
            {
                return {};
            }
        return rec_pol(gametes[g1], gametes[g2], mutations);
    }

    template <typename queue_type, typename mutation_model, typename mcont_t>
    std::vector<uint_t>
    generate_new_mutations(queue_type &recycling_bin, const gsl_rng *r,
                           const double &mu, mcont_t &mutations,
                           const std::size_t g, const mutation_model &mmodel)
    /// Return a vector of keys to new mutations.  The keys
    /// will be sorted according to mutation postition.
    {
        unsigned nm = gsl_ran_poisson(r, mu);
        std::vector<uint_t> rv;
        rv.reserve(nm);
        for (unsigned i = 0; i < nm; ++i)
            {
                rv.emplace_back(fwdpp_internal::mmodel_dispatcher(
                    mmodel, g, mutations, recycling_bin));
            }
        std::sort(rv.begin(), rv.end(),
                  [&mutations](const uint_t a, const uint_t b) {
                      return mutations[a].pos < mutations[b].pos;
                  });
        return rv;
    }

    template <typename gcont_t, typename mcont_t, typename queue_type>
    std::size_t
    mutate_recombine(
        const std::vector<uint_t> &new_mutations,
        const std::vector<double> &breakpoints,
        // const mutation_model &mmodel, const recombination_model &recmodel,
        const std::size_t g1, const std::size_t g2, gcont_t &gametes,
        mcont_t &mutations, queue_type &gamete_recycling_bin,
        typename gcont_t::value_type::mutation_container &neutral,
        typename gcont_t::value_type::mutation_container &selected)
    {
        if (new_mutations.empty() && breakpoints.empty())
            {
                return g1;
            }
        else if (breakpoints.empty())
            {
                neutral.clear();
                selected.clear();
                neutral.reserve(std::max(gametes[g1].mutations.size(),
                                         gametes[g2].mutations.size())
                                + new_mutations.size());
                selected.reserve(std::max(gametes[g1].smutations.size(),
                                          gametes[g2].smutations.size())
                                 + new_mutations.size());
                auto nb = gametes[g1].mutations.begin(),
                     sb = gametes[g1].smutations.begin();
                const auto ne = gametes[g1].mutations.end(),
                           se = gametes[g1].smutations.end();
                for (auto &&m : new_mutations)
                    {
                        if (mutations[m].neutral)
                            {
                                auto t = std::upper_bound(
                                    nb, ne, mutations[m].pos,
                                    [&mutations](
                                        const double &v,
                                        const std::size_t mut) noexcept {
                                        return v < mutations[mut].pos;
                                    });
                                neutral.insert(neutral.end(), nb, t);
                                neutral.push_back(m);
                                nb = t;
                            }
                        else
                            {
                                auto t = std::upper_bound(
                                    sb, se, mutations[m].pos,
                                    [&mutations](
                                        const double &v,
                                        const std::size_t mut) noexcept {
                                        return v < mutations[mut].pos;
                                    });
                                selected.insert(selected.end(), sb, t);
                                selected.push_back(m);
                                sb = t;
                            }
                    }
                neutral.insert(neutral.end(), nb, ne);
                selected.insert(selected.end(), sb, se);
// std::cerr << "returning a new mutant " <<
// new_mutations.size()
//         << '\n';
#ifndef NDEBUG
                std::unordered_map<uint_t, uint_t> n, s;
                for (auto &&ni : neutral)
                    n[ni]++;
                for (auto &&ni : selected)
                    s[ni]++;
                for (auto &&ni : n)
                    {
                        if (ni.second > 1)
                            {
                                std::cout << ni.second << ": -> ";
                                for (auto x : neutral)
                                    {
                                        std::cout << x << ' ';
                                    }
                                std::cout << '\n';
                                for (auto mi : new_mutations)
                                    {
                                        std::cout << mi << ' ';
                                    }
                                std::cout << '\n';
                            }
                        assert(ni.second == 1);
                    }
#endif
                return fwdpp_internal::recycle_gamete(
                    gametes, gamete_recycling_bin, neutral, selected);
            }
        neutral.clear();
        selected.clear();
        neutral.reserve(std::max(gametes[g2].mutations.size(),
                                 gametes[g2].mutations.size())
                        + new_mutations.size());
        selected.reserve(std::max(gametes[g2].smutations.size(),
                                  gametes[g2].smutations.size())
                         + new_mutations.size());

        assert(breakpoints.back() == std::numeric_limits<double>::max());
        assert(std::is_sorted(breakpoints.begin(), breakpoints.end()));
        assert(std::is_sorted(
            new_mutations.begin(), new_mutations.end(),
            [&mutations](const std::uint32_t i, const std::uint32_t j) {
                return mutations[i].pos < mutations[j].pos;
            }));

        auto itr = gametes[g1].mutations.cbegin();
        auto jtr = gametes[g2].mutations.cbegin();
        auto itr_s = gametes[g1].smutations.cbegin();
        auto jtr_s = gametes[g2].smutations.cbegin();
        auto itr_e = gametes[g1].mutations.cend();
        auto itr_s_e = gametes[g1].smutations.cend();
        auto jtr_e = gametes[g2].mutations.cend();
        auto jtr_s_e = gametes[g2].smutations.cend();

        auto next_mutation = new_mutations.cbegin();
        for (auto i = breakpoints.cbegin(); i != breakpoints.cend();)
            {
                if (next_mutation != new_mutations.cend()
                    && mutations[*next_mutation].pos < *i)
                    {
                        // next_pos = mutations[*next_mutation].pos;
                        // is_mut = true;
                        itr = fwdpp_internal::rec_gam_updater(
                            itr, itr_e, mutations, neutral,
                            mutations[*next_mutation].pos);
                        itr_s = fwdpp_internal::rec_gam_updater(
                            itr_s, itr_s_e, mutations, selected,
                            mutations[*next_mutation].pos);
                        jtr = fwdpp_internal::rec_update_itr(
                            jtr, jtr_e, mutations,
                            mutations[*next_mutation].pos);
                        jtr_s = fwdpp_internal::rec_update_itr(
                            jtr_s, jtr_s_e, mutations,
                            mutations[*next_mutation].pos);
                        if (mutations[*next_mutation].neutral)
                            {
                                neutral.push_back(*next_mutation);
                            }
                        else
                            {
                                selected.push_back(*next_mutation);
                            }
                        ++next_mutation;
                    }
                else
                    {
                        itr = fwdpp_internal::rec_gam_updater(
                            itr, itr_e, mutations, neutral, *i);
                        itr_s = fwdpp_internal::rec_gam_updater(
                            itr_s, itr_s_e, mutations, selected, *i);
                        jtr = fwdpp_internal::rec_update_itr(jtr, jtr_e,
                                                             mutations, *i);
                        jtr_s = fwdpp_internal::rec_update_itr(jtr_s, jtr_s_e,
                                                               mutations, *i);
                        std::swap(itr, jtr);
                        std::swap(itr_s, jtr_s);
                        std::swap(itr_e, jtr_e);
                        std::swap(itr_s_e, jtr_s_e);
                        ++i;
                    }
            }
        if (next_mutation != new_mutations.cend())
            {
                throw std::runtime_error("we didn't do all the mutations");
            }
#ifndef NDEBUG
        std::unordered_map<uint_t, uint_t> n, s;
        for (auto &&ni : neutral)
            n[ni]++;
        for (auto &&ni : selected)
            s[ni]++;
        for (auto &&ni : n)
            {
                if (ni.second > 1)
                    {
                        std::cout << ni.second << ": -> ";
                        for (auto x : neutral)
                            {
                                std::cout << x << ' ';
                            }
                        std::cout << '\n';
                        for (auto mi : new_mutations)
                            {
                                std::cout << mi << ' ';
                            }
                        std::cout << '\n';
                    }
                assert(ni.second == 1);
            }
        for (auto &&ni : s)
            {
                assert(ni.second == 1);
            }
#endif
        return fwdpp_internal::recycle_gamete(gametes, gamete_recycling_bin,
                                              neutral, selected);
    }
}

#endif
