#ifndef FWDPP_TS_GENERATE_OFFSPRING_HPP
#define FWDPP_TS_GENERATE_OFFSPRING_HPP

#include <tuple>
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

        template <typename genetic_param_holder, typename poptype,
                  typename recmodel, typename mutmodel, typename recycling_bin,
                  typename mutation_key_container>
        //TODO: replace with a struct?
        std::tuple<std::vector<double>, std::vector<uint_t>>
        generate_mutations_and_breakpoints(
            std::size_t parent, std::size_t parental_gamete,
            const recmodel& generate_breakpoints,
            const mutmodel& generate_mutations,
            recycling_bin& mutation_recycling_bin,
            recycling_bin& gamete_recycling_bin,
            mutation_key_container& neutral, mutation_key_container& selected,
            poptype& pop)
        {
            auto breakpoints = generate_breakpoints();
            auto new_mutation_keys = fwdpp_internal::mmodel_dispatcher(
                generate_mutations, pop.diploids[parent],
                pop.gametes[parental_gamete], pop.mutations,
                mutation_recycling_bin);
            return std::make_tuple(std::move(breakpoints),
                                   std::move(new_mutation_keys));
        }

        // NOTE The entire problem that we have here is due
        // to the idea that we may or may not want to put neutral mutations
        // into the gametes.  This function could be much simpler
        // if we just returned breakpoints and new mutation keys.
        //
        // We almost certainly DO want to reuse fwdpp::mutate_recombine,
        // which constrains us.  However, a lot of the hand-wringing is having
        // to do with not wanting to copy the vector of new mutation keys. But,
        // fwdpp::mutate_recombine uses an iterator range to update the new gamete
        // in one place, and a range-based for loop in another.  We should consider
        // a simplification where we change the API of that function to take a
        // pair of iterators. We could to likewise with breakpoints, too.
        template <typename genetic_param_holder,
                  typename mutation_handling_policy, typename poptype>
        inline TS_NODE_INT
        generate_offspring_details(
            fwdpp::sugar::SINGLELOC_TAG, const gsl_rng* r,
            const std::pair<std::size_t, std::size_t> parents,
            const ts_bookkeeper& bookkeeper,
            const mutation_handling_policy& mutation_policy, poptype& pop,
            genetic_param_holder& genetics,
            typename poptype::diploid_t& offspring, table_collection& tables)
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
            //TODO remove redundancy below
            //TODO the different order access into the tuple is confusing/error-prone
            auto pid = get_parent_ids(bookkeeper.first_parental_index,
                                      parents.first, swap1);
            auto gamete_data = generate_mutations_and_breakpoints(
                parents.first, p1g1, genetics.generate_breakpoints,
                genetics.generate_mutations, genetics.mutation_recycling_bin,
                genetics.gamete_recycling_bin, genetics.neutral,
                genetics.selected, pop);
            tables.add_offspring_data(
                bookkeeper.next_index, std::get<0>(gamete_data),
                std::get<1>(gamete_data), pid, bookkeeper.offspring_deme,
                bookkeeper.birth_time);
            process_new_mutations(std::get<1>(gamete_data), pop.mutations,
                                  mutation_policy);
            offspring.first = mutate_recombine(
                std::get<1>(gamete_data), std::get<0>(gamete_data), p1g1, p1g2,
                pop.gametes, pop.mutations, genetics.gamete_recycling_bin,
                genetics.neutral, genetics.selected);
            pid = get_parent_ids(bookkeeper.first_parental_index,
                                 parents.second, swap2);
            gamete_data = generate_mutations_and_breakpoints(
                parents.second, p2g1, genetics.generate_breakpoints,
                genetics.generate_mutations, genetics.mutation_recycling_bin,
                genetics.gamete_recycling_bin, genetics.neutral,
                genetics.selected, pop);
            tables.add_offspring_data(
                bookkeeper.next_index + 1, std::get<0>(gamete_data),
                std::get<1>(gamete_data), pid, bookkeeper.offspring_deme,
                bookkeeper.birth_time);
            process_new_mutations(std::get<1>(gamete_data), pop.mutations,
                                  mutation_policy);
            offspring.first = mutate_recombine(
                std::get<1>(gamete_data), std::get<0>(gamete_data), p2g1, p2g2,
                pop.gametes, pop.mutations, genetics.gamete_recycling_bin,
                genetics.neutral, genetics.selected);
            return bookkeeper.next_index + 2;
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
        {
            return generate_offspring_details(
                typename poptype::popmodel_t(), r, parents, bookkeeper,
                mutation_policy, pop, genetics, offspring, tables);
        }
    } // namespace ts
} // namespace fwdpp

#endif
