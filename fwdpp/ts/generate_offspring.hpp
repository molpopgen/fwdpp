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

        template <typename genetic_param_holder,
                  typename mutation_handling_policy, typename poptype>
        inline TS_NODE_INT
        generate_offspring(const gsl_rng* r,
                           const std::pair<std::size_t, std::size_t> parents,
                           const ts_bookkeeper& bookkeeper,
                           const mutation_handling_policy& mutation_policy,
                           poptype& pop, genetic_param_holder& genetics,
                           typename poptype::diploid_t& offspring,
                           table_collection& tables)
		//TODO: dispatch this out differently for slocuspop vs mlocuspop
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
			//TODO remove redundancy below
            auto p1id = get_parent_ids(bookkeeper.first_parental_index,
                                       parents.first, swap1);
            auto breakpoints = genetics.generate_breakpoints();
            auto new_mutation_keys = fwdpp_internal::mmodel_dispatcher(
                genetics.generate_mutations, pop.diploids[parents.first],
                pop.gametes[p1g1], pop.mutations,
                genetics.mutation_recycling_bin);
            tables.add_offspring_data(
                bookkeeper.next_index, breakpoints, new_mutation_keys, p1id,
                bookkeeper.offspring_deme, bookkeeper.birth_time);
            process_new_mutations(new_mutation_keys, pop.mutations,
                                  mutation_policy);
            offspring.first = mutate_recombine(
                new_mutation_keys, breakpoints, p1g1, p1g2, pop.gametes,
                pop.mutations, genetics.gamete_recycling_bin, genetics.neutral,
                genetics.selected);
            auto p2id = get_parent_ids(bookkeeper.first_parental_index,
                                       parents.second, swap2);
            breakpoints = genetics.generate_breakpoints();
            new_mutation_keys = fwdpp_internal::mmodel_dispatcher(
                genetics.generate_mutations, pop.diploids[parents.second],
                pop.gametes[p2g1], pop.mutations,
                genetics.mutation_recycling_bin);
            tables.add_offspring_data(
                bookkeeper.next_index + 1, breakpoints, new_mutation_keys,
                p2id, bookkeeper.offspring_deme, bookkeeper.birth_time);
            process_new_mutations(new_mutation_keys, pop.mutations,
                                  mutation_policy);
            offspring.first = mutate_recombine(
                new_mutation_keys, breakpoints, p2g1, p2g2, pop.gametes,
                pop.mutations, genetics.gamete_recycling_bin, genetics.neutral,
                genetics.selected);
            return bookkeeper.next_index + 2;
        }
    } // namespace ts
} // namespace fwdpp

#endif
