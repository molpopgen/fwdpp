#ifndef FWDPP_TS_RECYCLING_HPP
#define FWDPP_TS_RECYCLING_HPP

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <vector>
#include <limits>
#include <cstdint>
#include <deque>
#include <queue>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/simfunctions/recycling.hpp>

namespace fwdpp
{
    namespace ts
    {
        inline flagged_mutation_queue
        make_mut_queue(
            const std::vector<std::uint32_t> &mcounts,
            const std::vector<std::uint32_t> &counts_from_preserved_nodes)
        /// \brief Make a mutation recycling queue for simulations with tree sequences
        /// \param mcounts Contribution of extant nodes to mutation counts
        /// \param counts_from_preserved_nodes Contribution of extinct nodes to mutation counts
        ///
        /// \version 0.7.0 Added to fwdpp
        ///
        /// \returns flagged_mutation_queue
        {
            flagged_mutation_queue::value_type mutation_recycling_bin;
            for (std::size_t i = 0; i < mcounts.size(); ++i)
                {
                    if (mcounts[i] + counts_from_preserved_nodes[i] == 0)
                        {
                            mutation_recycling_bin.push(i);
                        }
                }
            return flagged_mutation_queue(std::move(mutation_recycling_bin));
        }

        inline flagged_mutation_queue
        make_mut_queue(
            const std::vector<std::size_t> &preserved_mutation_indexes,
            const std::size_t num_mutations)
        /// \brief Make a mutation recycling queue for simulations with tree sequences
        /// \param preserved_mutation_indexes Vector of preserved mutation indexes returned by simplification
        /// \param num_mutations The total number of mutations currently allocated in the population
        ///
        /// \returns flagged_mutation_queue
        ///
        /// \version 0.7.3 Added to fwdpp
        ///
        /// This overload may be preferable to the other when the following conditions apply:
        /// 1. Mutation counts in the entire simulation are not of interest during the simulation.
        /// 2. There is no need/wish to remove fixations from the gametes/tables during the simulation.
        /// 3. There are large numbers of ancient samples being recorded.
        ///
        /// The first two conditions are required for correct results.  The third condition is optional,
        /// but big speedups will be seen for that case.
        ///
        /// This function generates a recycling queue using a fast O(N) algorithm.
        /// This, it will outperform the other overload based
        /// on tree traversal when the number of trees is very large, as is the case when large numbers
        /// of ancestral samples are registered during a simulation.
        ///
        {
            std::vector<std::size_t> mindexes(num_mutations);
            std::iota(mindexes.begin(), mindexes.end(), 0);
            for (auto i : preserved_mutation_indexes)
                {
                    mindexes[i] = std::numeric_limits<std::size_t>::max();
                }
            flagged_mutation_queue::value_type rv;
            for (auto i : mindexes)
                {
                    if (i != std::numeric_limits<std::size_t>::max())
                        {
                            rv.push(i);
                        }
                }
            return flagged_mutation_queue(std::move(rv));
        }

        namespace detail
        {
            template <typename mcont_t, typename lookup_table>
            void
            process_mutation_index(mcont_t &mutations, lookup_table &lookup,
                                   const std::size_t i)
            /// Implementation detail for flag_mutations_for_recycling
            {
                auto itr = lookup.equal_range(mutations[i].pos);
                mutations[i].pos = std::numeric_limits<double>::max();
                while (itr.first != itr.second)
                    {
                        if (itr.first->second == i)
                            {
                                lookup.erase(itr.first);
                                break;
                            }
                        ++itr.first;
                    }
            }

            template <typename mcont_t, typename mutation_count_container,
                      typename lookup_table>
            inline void
            process_fixations(mcont_t &mutations,
                              mutation_count_container &mcounts,
                              mcont_t & /*fixations*/,
                              std::vector<uint_t> & /*fixation_times*/,
                              lookup_table &mut_lookup,
                              const uint_t /*generation*/, const std::size_t i,
                              std::false_type, std::false_type)
            // Remove all fixations.
            // Do not record fixations
            {
                process_mutation_index(mutations, mut_lookup, i);
                mcounts[i] = 0;
            }

            template <typename mcont_t, typename mutation_count_container,
                      typename lookup_table>
            inline void
            process_fixations(mcont_t &mutations,
                              mutation_count_container &mcounts,
                              mcont_t &fixations,
                              std::vector<uint_t> &fixation_times,
                              lookup_table &mut_lookup,
                              const uint_t generation, const std::size_t i,
                              std::false_type, std::true_type)
            // Remove all fixations
            // Record fixations
            {
                fixations.push_back(mutations[i]);
                fixation_times.push_back(generation);
                process_mutation_index(mutations, mut_lookup, i);
                mcounts[i] = 0;
            }

            struct equal_range_comparison
            {
                template <typename mutation_type>
                inline bool
                operator()(const mutation_type &m, const double p) const
                {
                    return m.pos < p;
                }
                template <typename mutation_type>
                inline bool
                operator()(const double p, const mutation_type &m) const
                {
                    return p < m.pos;
                }
            };

            template <typename mcont_t, typename mutation_count_container,
                      typename lookup_table>
            inline void
            process_fixations(mcont_t &mutations,
                              mutation_count_container &mcounts,
                              mcont_t & /*fixations*/,
                              std::vector<uint_t> & /*fixation_times*/,
                              lookup_table &mut_lookup,
                              const uint_t /*generation*/, const std::size_t i,
                              std::true_type, std::false_type)
            // Only remove neutral fixations
            // Do not record fixations
            {
                if (mutations[i].neutral)
                    {
                        process_mutation_index(mutations, mut_lookup, i);
                        mcounts[i] = 0;
                    }
            }

            template <typename mcont_t, typename mutation_count_container,
                      typename lookup_table>
            inline void
            process_fixations(mcont_t &mutations,
                              mutation_count_container &mcounts,
                              mcont_t &fixations,
                              std::vector<uint_t> &fixation_times,
                              lookup_table &mut_lookup,
                              const uint_t generation, const std::size_t i,
                              std::true_type, std::true_type)
            // Only remove neutral fixations.
            // Record fixations.
            // The recording is in order of mutation position.
            {
                if (mutations[i].neutral)
                    {
                        // Record the fixation here,
                        // as we will recycle this variant
                        auto loc = std::lower_bound(
                            fixations.begin(), fixations.end(),
                            mutations[i].pos,
                            [](const typename mcont_t::value_type &m,
                               const double val) { return m.pos < val; });
                        auto d = std::distance(fixations.begin(), loc);
                        fixations.insert(loc, mutations[i]);
                        fixation_times.insert(fixation_times.begin() + d,
                                              generation);
                        process_mutation_index(mutations, mut_lookup, i);
                        mcounts[i] = 0;
                    }
                else
                    {
                        // The logic here is tougher.  We are preserving
                        // selected mutations, and thus we only
                        // want to record their fixations once. Thus,
                        // we have to make sure that this fixation
                        // hasn't been recorded at all.
                        auto loc_range = std::equal_range(
                            fixations.begin(), fixations.end(),
                            mutations[i].pos,
                            detail::equal_range_comparison());
                        if (std::find(loc_range.first, loc_range.second,
                                      mutations[i])
                            == loc_range.second)
                            {
                                auto d = std::distance(fixations.begin(),
                                                       loc_range.first);
                                fixations.insert(loc_range.first,
                                                 mutations[i]);
                                fixation_times.insert(
                                    fixation_times.begin() + d, generation);
                            }
                    }
            }
        } // namespace detail

        template <typename poptype, typename mutation_count_container,
                  typename preserve_selected_fixations,
                  typename record_fixations>
        void
        flag_mutations_for_recycling(
            poptype &pop,
            mutation_count_container &mcounts_from_preserved_nodes,
            const uint_t twoN, const uint_t generation,
            const preserve_selected_fixations preserve,
            const record_fixations record)
        /*! Mark mutations for recycling.
         *
         * This function should be called immediately after simplification and mutation counting.
         *
         * \param pop A mutation container
         * \param mcounts_from_preserved_nodes A container recording the counts of each mutation in ancient samples
         * \param twoN Twice the current population size
         * \param generation Current generation/time step in the simulation
         * \param preserve Policy for handling selected fixations.  See below.
         * \param record Policy for recording fixation events. See below.
         *
         * A mutation is marked for recycling if one of the following conditions holds:
         * 1. The sum of \a pop.mcounts and \a mcounts_from_preserved_nodes is zero.
         * 2. \a pop.mcounts == \a twoN, mcounts_from_preserved_nodes is zero, 
         * and \a record is std::false_type or the mutation is neutral.
         *
         * Condition 1 refers to extinct mutations and condition 2 refers to fixations.
         *
         * When \a preserve is std::true_type, selected fixations are retained
         * in the population.  We do this because simulations of phenotypes (as opposed
         * to relative fitness) require tracking the contribution of fixation to 
         * trait values.
         *
         * All variants matching the above criteria have their record in \a pop.mcounts
         * set to zero and their position is removed from \a pop.mut_lookup.  Further, the 
         * mutation's position in \a pop.mutations is set to std::numeric_limits<double>::max(),
         * to signify an "invalid" mutation.
         *
         * If \a record is std::false_type, no fixation recording takes place.  However, if
         * \a record is std::true_type, a record is entered into \a pop.fixations and
         * \a pop.fixation_times.  If \a preserve is also std::true_type, the recording is 
         * slightly more expensive because we have to guard against repeated recording.  For this
         * case, the fixations and fixation times containers in \a pop are kept sorted by fixation position.
         *
         * An advanced use of this function is to create branchless code with respect to 
         * fixation handling.  One can generate a closure wrapping this function with the desired
         * values of \a preserve and \a record at the start of a simulation, thus avoiding a bunch of 
         * "if" statements prior to each call.
         *
         * \note There are several caveats to fixation recording (\a record == std::true_type).  First,
         * it is only as accurate as the time steps between simplifications is small.  \a pop will only 
         * contain exact fixation times if simplification occurs every generation.  The second limitation
         * is that it will not give correct results if ancient samples are being tracked in a simulation.
         * This has to do with a "todo" item for this function described below.  Thus, fixation recording 
         * with this function should be considered limited to the case of regular simplification and no
         * ancient sample preservation.
         *
         * \warning The \a preserve argument needs to be consistent with the desired behavior of
         * ts::remove_fixations_from_gametes.
         *
         * \version 0.7.0 Added to library
         * \version 0.7.1 Updated to change recycled mutation positions to max value of a double.
         * Refactor API to take compile-time policies for fixation handling.
         *
         * \todo Improve treatment of fixations by allowing for variants fixed in alive AND 
         * ancient samples to be flagged.
         * \todo Return a recycling queue?
         */
        {
            for (std::size_t i = 0; i < pop.mcounts.size(); ++i)
                {
                    if (mcounts_from_preserved_nodes[i] == 0)
                        {
                            if (pop.mcounts[i] > twoN)
                                {
                                    throw std::runtime_error(
                                        "mutation count out of range");
                                }
                            if (pop.mcounts[i] == twoN)
                                {
                                    detail::process_fixations(
                                        pop.mutations, pop.mcounts,
                                        pop.fixations, pop.fixation_times,
                                        pop.mut_lookup, generation, i,
                                        preserve, record);
                                }
                            else if (pop.mcounts[i] == 0)
                                {
                                    detail::process_mutation_index(
                                        pop.mutations, pop.mut_lookup, i);
                                }
                        }
                }
        }
    } // namespace ts
} // namespace fwdpp

#endif
