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
#include <fwdpp/internal/mutation_internal.hpp>

namespace fwdpp
{
    namespace ts
    {
        struct mut_rec_intermediates
        /// Contains the output of recombination and mutation when generating offspring gametes.
        /// Objects of this type are returned by fwdpp::ts::generate_offspring.
        /// \version 0.7.4 Added to fwdpp::ts.
        {
            /// 1 if the parent gametes were swapped, 0 otherwise
            int swapped;
            /// Recombination breakpoints
            std::vector<double> breakpoints;
            /// Keys to new mutations
            std::vector<uint_t> mutation_keys;

            template <typename B, typename M>
            mut_rec_intermediates(int s, B&& b, M&& m)
                : swapped{ s }, breakpoints{ std::forward<B>(b) },
                  mutation_keys{ std::forward<M>(m) }
            {
            }
        };

        struct all_mutations
        /// Dispatch tag. See docs for fwdpp::ts::generate_offspring
        /// \version 0.7.4 Added to fwdpp::ts.
        {
        };

        struct selected_variants_only
        /// Dispatch tag. See docs for fwdpp::ts::generate_offspring
        /// \version 0.7.4 Added to fwdpp::ts.
        {
        };

        namespace detail
        {

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
                const parental_data parent,
                const recmodel& generate_breakpoints,
                const mutmodel& generate_mutations,
                const mutation_handling_policy& mutation_policy,
                recycling_bin& mutation_recycling_bin,
                recycling_bin& gamete_recycling_bin,
                mutation_key_container& neutral,
                mutation_key_container& selected, poptype& pop)
            {
                auto gamete_data = generate_mutations_and_breakpoints(
                    parent.index, parent.gamete1, parent.swapped,
                    generate_breakpoints, generate_mutations,
                    mutation_recycling_bin, pop);
                auto range = process_new_mutations(
                    gamete_data.mutation_keys, pop.mutations, mutation_policy);
                std::size_t offspring_gamete = mutate_recombine(
                    range, gamete_data.breakpoints, parent.gamete1,
                    parent.gamete2, pop.gametes, pop.mutations,
                    gamete_recycling_bin, neutral, selected);
                return std::make_pair(offspring_gamete,
                                      std::move(gamete_data));
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
            generate_offspring_details(
                fwdpp::poptypes::MULTILOC_TAG, const gsl_rng* r,
                const std::pair<std::size_t, std::size_t> parents,
                const mutation_handling_policy& mutation_policy, poptype& pop,
                genetic_param_holder& genetics,
                typename poptype::diploid_t& offspring)
            {
                int swap1 = (gsl_rng_uniform(r) < 0.5) ? 1 : 0;
                int swap2 = (gsl_rng_uniform(r) < 0.5) ? 1 : 0;
                int ttl_swaps_1 = swap1;
                int ttl_swaps_2 = swap2;

                decltype(mut_rec_intermediates::breakpoints) all_breakpoints_1,
                    all_breakpoints_2;
                decltype(mut_rec_intermediates::mutation_keys) all_mut_keys_1,
                    all_mut_keys_2;

                const auto& irec = genetics.interlocus_recombination;
                for (std::size_t i = 0; i < offspring.size; ++i)
                    {
                        if (i > 0)
                            {
                                // between-locus rec, parent 1
                                ttl_swaps_1 += irec[i - 1]();
                                // between-locus rec, parent 2
                                ttl_swaps_2 += irec[i - 1]();
                            }
                        auto p1g1 = pop.diploids[parents.first][i].first;
                        auto p1g2 = pop.diploids[parents.first][i].second;
                        auto p2g1 = pop.diploids[parents.second][i].first;
                        auto p2g2 = pop.diploids[parents.second][i].second;
                        if (ttl_swaps_1 % 2 != 0.)
                            {
                                std::swap(p1g1, p1g2);
                                if (i > 0)
                                    {
                                        all_breakpoints_1.push_back(
                                            pop.locus_boundaries[i - 1]
                                                .second);
                                    }
                            }
                        if (ttl_swaps_2 % 2 != 0.)
                            {
                                std::swap(p2g1, p2g2);
                                if (i > 0)
                                    {
                                        all_breakpoints_2.push_back(
                                            pop.locus_boundaries[i - 1]
                                                .second);
                                    }
                            }
                        auto first_gamete_data = generate_offspring_gamete(
                            parental_data{ parents.first, p1g1, p1g2, swap1 },
                            genetics.generate_breakpoints,
                            genetics.generate_mutations, mutation_policy,
                            genetics.mutation_recycling_bin,
                            genetics.gamete_recycling_bin, genetics.neutral,
                            genetics.selected, pop);
                        auto second_gamete_data = generate_offspring_gamete(
                            parental_data{ parents.second, p2g1, p2g2, swap2 },
                            genetics.generate_breakpoints,
                            genetics.generate_mutations, mutation_policy,
                            genetics.mutation_recycling_bin,
                            genetics.gamete_recycling_bin, genetics.neutral,
                            genetics.selected, pop);
                        // Update the offspring's gametes.
                        offspring[i].first = first_gamete_data.first;
                        offspring[i].second = second_gamete_data.first;
                        // Add mutations to the return values
                        all_mut_keys_1.insert(
                            end(all_mut_keys_1),
                            begin(first_gamete_data.mutation_keys),
                            end(first_gamete_data.mutation_keys));
                        all_mut_keys_2.insert(
                            end(all_mut_keys_2),
                            begin(second_gamete_data.mutation_keys),
                            end(second_gamete_data.mutation_keys));
                        // Update recombination breakpoints
                        if (!first_gamete_data.breakpoints.empty())
                            {
                                ttl_swaps_1
                                    += first_gamete_data.breakpoints.size()
                                       - 1;
                                all_breakpoints_1.insert(
                                    end(all_breakpoints_1),
                                    begin(first_gamete_data.breakpoints),
                                    end(first_gamete_data.breakpoints) - 1);
                            }
                        if (!second_gamete_data.breakpoints.empty())
                            {
                                ttl_swaps_2
                                    += second_gamete_data.breakpoints.size()
                                       - 1;
                                all_breakpoints_2.insert(
                                    end(all_breakpoints_1),
                                    begin(second_gamete_data.breakpoints),
                                    end(second_gamete_data.breakpoints) - 1);
                            }
                    }

                return std::make_pair(
                    mut_rec_intermediates(swap1, std::move(all_breakpoints_1),
                                          std::move(all_breakpoints_2)),
                    mut_rec_intermediates(swap2, std::move(all_breakpoints_2),
                                          std::move(all_mut_keys_2)));
            }

        } // namespace detail

        template <typename genetic_param_holder,
                  typename mutation_handling_policy, typename poptype>
        std::pair<mut_rec_intermediates, mut_rec_intermediates>
        generate_offspring(const gsl_rng* r,
                           const std::pair<std::size_t, std::size_t> parents,
                           const mutation_handling_policy& mutation_policy,
                           poptype& pop, genetic_param_holder& genetics,
                           typename poptype::diploid_t& offspring)
        /// \brief Generate offspring gametes and return breakpoints plus mutation keys
        ///
        /// \param r Random number generator
        /// \param parents Indexes of the offspring parents in \a pop
        /// \param mutation_policy Either all_mutations or selected_variants_only.  See below.
        /// \param pop Either fwdpp::poptypes::slocuspop or fwdpp::poptypes::mlocuspop.
        /// \param genetics A duck type of fwdpp::genetic_parameters.
        /// \param offspring The offspring for which we will generate gametes.
        ///
        /// The final three parameters will be modified.
        ///
        /// \return A pair of mut_rec_intermediates, corresponding to what is passed on from
        /// each parent.
        ///
        /// The operations are dispatched out to functions based on the population type.  These
        /// functions make calls to fwdpp::mutate_recombine.
        ///
        /// The parameter \a mutation_policy governs what types of mutations are entered
        /// into the gametes of \a offspring.  If the policy is fwdpp::ts::all_mutations, then
        /// neutral and selected variants are placed in an offspring's gametes.  If the policy is
        /// fwdpp::ts::selected_variants_only, then only selected mutations are placed into the gametes.
        /// Regardless of the policy, ALL mutations are contained in the return value, with the idea
        /// that the caller will record them into a fwdpp::ts::table_collection.
        ///
        /// \note The function type genetic_param_holder::generate_mutations must return
        /// std::vector<fwdpp::uint_t>, with the values reflecting the locations of new
        /// mutations on \a pop.mutations.  Further, this function must already bind
        /// any relevant mutation rates (as this generate_offspring does not accept them
        /// as arguments).
        ///
        /// \version 0.7.4 Added to fwdpp::ts.
        ///
        {
            return detail::generate_offspring_details(
                typename poptype::popmodel_t(), r, parents, mutation_policy,
                pop, genetics, offspring);
        }
    } // namespace ts
} // namespace fwdpp

#endif
