#ifndef FWDPP_TS_GENERATE_OFFSPRING_HPP
#define FWDPP_TS_GENERATE_OFFSPRING_HPP

#include <tuple>
#include <vector>
#include <utility>
#include <type_traits>
#include <gsl/gsl_rng.h>
#include <fwdpp/wrapped_range.hpp>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/mutate_recombine.hpp>
#include <fwdpp/ts/get_parent_ids.hpp>
#include <fwdpp/ts/definitions.hpp>
#include <fwdpp/ts/table_collection.hpp>
#include <fwdpp/internal/mutation_internal.hpp>

namespace fwdpp
{
    namespace ts
    {
        struct mut_rec_intermediates
        {
            int swapped;
            std::vector<double> breakpoints;
            std::vector<uint_t> mutation_keys;
            template <typename B, typename M>
            mut_rec_intermediates(int s, B&& b, M&& m)
                : swapped{ s }, breakpoints{ std::forward<B>(b) },
                  mutation_keys{ std::forward<M>(m) }
            {
            }
        };

        struct all_mutations
        {
        };

        struct selected_variants_only
        {
        };

        //TODO: much of the stuff below can go into an internal namespace
        template <typename key_vector, typename mcont_t>
        inline wrapped_range<typename key_vector::iterator>
        process_new_mutations(const key_vector& new_mutation_keys,
                              const mcont_t&, all_mutations)
        {
            return make_wrapped_range(begin(new_mutation_keys),
                                      end(new_mutation_keys));
        }

        template <typename key_vector, typename mcont_t>
        inline wrapped_range<typename key_vector::iterator>
        process_new_mutations(key_vector& new_mutation_keys,
                              mcont_t& mutations, selected_variants_only)
        {
            auto itr = std::stable_partition(
                begin(new_mutation_keys), end(new_mutation_keys),
                [&mutations](typename key_vector::value_type k) {
                    return mutations[k].neutral == false;
                });
            return make_wrapped_range(begin(new_mutation_keys), itr);
        }

        template <typename poptype, typename recmodel, typename mutmodel,
                  typename recycling_bin>
        inline mut_rec_intermediates
        generate_mutations_and_breakpoints(
            std::size_t parent, std::size_t parental_gamete, int swapped,
            const recmodel& generate_breakpoints,
            const mutmodel& generate_mutations,
            recycling_bin& mutation_recycling_bin, poptype& pop)
        {
            auto breakpoints = generate_breakpoints();
            auto new_mutation_keys = fwdpp_internal::mmodel_dispatcher(
                generate_mutations, pop.diploids[parent],
                pop.gametes[parental_gamete], pop.mutations,
                mutation_recycling_bin);
            return mut_rec_intermediates(swapped, std::move(breakpoints),
                                         std::move(new_mutation_keys));
        }

        struct parental_data
        {
            std::size_t index, gamete1, gamete2;
            int swapped;
        };

        template <typename poptype, typename recmodel, typename mutmodel,
                  typename recycling_bin, typename mutation_key_container,
                  typename mutation_handling_policy>
        inline std::pair<std::size_t, mut_rec_intermediates>
        generate_offspring_gamete(
            const parental_data parent, const recmodel& generate_breakpoints,
            const mutmodel& generate_mutations,
            const mutation_handling_policy& mutation_policy,
            recycling_bin& mutation_recycling_bin,
            recycling_bin& gamete_recycling_bin,
            mutation_key_container& neutral, mutation_key_container& selected,
            poptype& pop)
        {
            auto gamete_data = generate_mutations_and_breakpoints(
                parent.index, parent.gamete1, parent.swapped,
                generate_breakpoints, generate_mutations,
                mutation_recycling_bin, pop);
            auto range = process_new_mutations(gamete_data.mutation_keys,
                                               pop.mutations, mutation_policy);
            std::size_t offspring_gamete = mutate_recombine(
                range, gamete_data.breakpoints, parent.gamete1, parent.gamete2,
                pop.gametes, pop.mutations, gamete_recycling_bin, neutral,
                selected);
            return std::make_pair(offspring_gamete, std::move(gamete_data));
        }

        template <typename genetic_param_holder,
                  typename mutation_handling_policy, typename poptype>
        std::pair<mut_rec_intermediates, mut_rec_intermediates>
        generate_offspring_details(
            fwdpp::poptypes::SINGLELOC_TAG, const gsl_rng* r,
            const std::pair<std::size_t, std::size_t> parents,
            const mutation_handling_policy& mutation_policy, poptype& pop,
            genetic_param_holder& genetics,
            typename poptype::diploid_t& offspring)
        {
            static_assert(
                std::is_same<typename std::remove_const<decltype(
                                 genetics.interlocus_recombination)>::type,
                             std::nullptr_t>::value,
                "invalid genetics policies detected");
            auto p1g1 = pop.diploids[parents.first].first;
            auto p1g2 = pop.diploids[parents.first].second;
            auto p2g1 = pop.diploids[parents.second].first;
            auto p2g2 = pop.diploids[parents.second].second;

            int swap1 = (gsl_rng_uniform(r) < 0.5) ? 1 : 0;
            int swap2 = (gsl_rng_uniform(r) < 0.5) ? 1 : 0;

            if (swap1)
                {
                    std::swap(p1g1, p1g2);
                }
            if (swap2)
                {
                    std::swap(p2g1, p2g2);
                }
            auto offspring_first_gamete_data = generate_offspring_gamete(
                parental_data{ parents.first, p1g1, p1g2, swap1 },
                genetics.generate_breakpoints, genetics.generate_mutations,
                mutation_policy, genetics.mutation_recycling_bin,
                genetics.gamete_recycling_bin, genetics.neutral,
                genetics.selected, pop);
            auto offspring_second_gamete_data = generate_offspring_gamete(
                parental_data{ parents.second, p2g1, p2g2, swap2 },
                genetics.generate_breakpoints, genetics.generate_mutations,
                mutation_policy, genetics.mutation_recycling_bin,
                genetics.gamete_recycling_bin, genetics.neutral,
                genetics.selected, pop);
            // Update the offspring's gametes.
            offspring.first = offspring_first_gamete_data.first;
            offspring.second = offspring_second_gamete_data.first;
            return std::make_pair(
                std::move(offspring_first_gamete_data.second),
                std::move(offspring_second_gamete_data.second));
        }

        template <typename genetic_param_holder,
                  typename mutation_handling_policy, typename poptype>
        std::pair<mut_rec_intermediates, mut_rec_intermediates>
        generate_offspring(const gsl_rng* r,
                           const std::pair<std::size_t, std::size_t> parents,
                           const mutation_handling_policy& mutation_policy,
                           poptype& pop, genetic_param_holder& genetics,
                           typename poptype::diploid_t& offspring)
        {
            return generate_offspring_details(typename poptype::popmodel_t(),
                                              r, parents, mutation_policy, pop,
                                              genetics, offspring);
        }
    } // namespace ts
} // namespace fwdpp

#endif