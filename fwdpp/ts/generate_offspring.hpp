#ifndef FWDPP_TS_GENERATE_OFFSPRING_HPP
#define FWDPP_TS_GENERATE_OFFSPRING_HPP

#include <type_traits>
#include <gsl/gsl_rng.h>
#include <fwdpp/mutate_recombine.hpp>
#include <fwdpp/ts/get_parent_ids.hpp>
#include <fwdpp/ts/definitions.hpp>
#include <fwdpp/ts/table_collection.hpp>
#include <fwdpp/internal/mutation_internal.hpp>

namespace fwdpp
{
    namespace ts
    {
        struct ts_bookkeeper
        {
            double birth_time;
            std::int32_t offspring_deme;
            TS_NODE_INT first_parental_index, next_index;
        };

        struct all_mutations
        {
        };
        struct selected_variants_only
        {
        };

        template <typename key_vector, typename mcont_t>
        inline void
        process_new_mutations(key_vector&, mcont_t&, all_mutations)
        {
        }

        template <typename key_vector, typename mcont_t>
        inline void
        process_new_mutations(key_vector& new_mutation_keys,
                              mcont_t& mutations, selected_variants_only)
        {
            auto itr = std::remove_if(
                new_mutation_keys.begin(), new_mutation_keys.end(),
                [&mutations](typename key_vector::value_type k) {
                    return mutations[k].neutral == true;
                });
            new_mutation_keys.erase(itr, new_mutation_keys.end());
        }

        template <typename breakpoint_function, typename new_mutation_fuction,
                  typename mutation_handling_policy, typename poptype,
                  typename mrecbin, typename grecbin>
        inline TS_NODE_INT
        generate_offspring(const gsl_rng* r,
                           const std::pair<std::size_t, std::size_t> parents,
                           const ts_bookkeeper& bookkeeper,
                           const breakpoint_function& recombination_function,
                           const new_mutation_fuction& mutation_function,
                           const mutation_handling_policy& mutation_policy,
                           poptype& pop,
                           typename poptype::diploid_t& offspring,
                           table_collection& tables,
                           mrecbin& mutation_recycling_bin,
                           grecbin& gamete_recycling_bin)
        {
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
            auto p1id = get_parent_ids(bookkeeper.first_parental_index,
                                       parents.first, swap1);
            auto p2id = get_parent_ids(bookkeeper.first_parental_index,
                                       parents.second, swap2);
            auto breakpoints = recombination_function();
            auto new_mutation_keys = fwdpp_internal::mmodel_dispatcher(
                mutation_function, pop.diploids[parents.first],
                pop.gametes[pop.diploids[parents.first].first], pop.mutations,
                mutation_recycling_bin);
            tables.add_offspring_data(
                bookkeeper.next_index, breakpoints, new_mutation_keys, p1id,
                bookkeeper.offspring_deme, bookkeeper.birth_time);
            process_new_mutations(new_mutation_keys, pop.mutations,
                                  mutation_policy);
            offspring.first = mutate_recombine(
                new_mutation_keys, breakpoints, p1g1, p1g2, pop.gametes,
                pop.mutations, gamete_recycling_bin, pop.neutral,
                pop.selected);
        }
    } // namespace ts
} // namespace fwdpp

#endif
