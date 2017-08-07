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

    template <typename mutation_queue_type, typename gamete_queue_type,
              typename gcont_t, typename mcont_t, typename mutation_model,
              typename recombination_model>
    std::size_t
    mutate_recombine(
        const gsl_rng *r, const double mutation_rate,
        const mutation_model &mmodel, const recombination_model &recmodel,
        const std::size_t g1, const std::size_t g2, gcont_t &gametes,
        mcont_t &mutations, mutation_queue_type &mut_recycling_bin,
        gamete_queue_type &gamete_recycling_bin,
        typename gcont_t::value_type::mutation_container &neutral,
        typename gcont_t::value_type::mutation_container &selected)
    {
        // The order here is to keep results same as previous fwdpp versions,
        // so we generate recombination breakpoints first, and then the new
        // mutations:
        auto breakpoints = recmodel(gametes[g1], gametes[g2], mutations);
        auto new_mutations = generate_new_mutations(
            mut_recycling_bin, r, mutation_rate, mutations, g1, mmodel);

        if (new_mutations.empty() && breakpoints.empty())
            {
                return g1;
            }
        else if (breakpoints.empty())
            {
                neutral.insert(neutral.end(), gametes[g1].mutations.begin(),
                               gametes[g1].mutations.end());
                selected.insert(selected.end(), gametes[g1].smutations.begin(),
                                gametes[g1].smutations.end());
                for (auto &&m : new_mutations)
                    {
                        auto c = (mutations[m].neutral) ? &neutral : &selected;
                        c->emplace(std::upper_bound(
                                       c->begin(), c->end(), mutations[m].pos,
                                       [&mutations](
                                           const double &__value,
                                           const std::size_t __mut) noexcept {
                                           assert(__mut < mutations.size());
                                           return __value
                                                  < mutations[__mut].pos;
                                       }),
                                   m);
                    }
                std::cerr << "returning a new mutant " << new_mutations.size()
                          << '\n';
                return fwdpp_internal::recycle_gamete(
                    gametes, gamete_recycling_bin, neutral, selected);
            }
        std::cerr << breakpoints.size() << ' ' << new_mutations.size() << '\n';

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
                // Get the next relevant position
                bool is_mut = false;
                double next_pos = *i;
                if (next_mutation != new_mutations.cend()
                    && mutations[*next_mutation].pos < *i)
                    {
                        next_pos = mutations[*next_mutation].pos;
                        is_mut = true;
                    }
                std::cerr << is_mut << ' ' << next_pos << ' '
                          << std::distance(next_mutation, new_mutations.cend())
                          << ' ' << std::distance(i, breakpoints.cend())
                          << '\n';
                itr = fwdpp_internal::rec_gam_updater(itr, itr_e, mutations,
                                                      neutral, next_pos);
                itr_s = fwdpp_internal::rec_gam_updater(
                    itr_s, itr_s_e, mutations, selected, next_pos);
                jtr = fwdpp_internal::rec_update_itr(jtr, jtr_e, mutations,
                                                     next_pos);
                jtr_s = fwdpp_internal::rec_update_itr(jtr_s, jtr_s_e,
                                                       mutations, next_pos);

                std::swap(itr, jtr);
                std::swap(itr_s, jtr_s);
                std::swap(itr_e, jtr_e);
                std::swap(itr_s_e, jtr_s_e);
                if (is_mut)
                    {
                        std::cerr << "is a mutation!\n";
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
                        ++i;
                    }
            }
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
        std::cerr << "DONE\n";
        return fwdpp_internal::recycle_gamete(gametes, gamete_recycling_bin,
                                              neutral, selected);
    }
}

#endif
