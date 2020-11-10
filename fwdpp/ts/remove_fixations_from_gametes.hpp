#ifndef FWDPP_TS_REMOVE_FIXATIONS_FROM_GAMETES_HPP
#define FWDPP_TS_REMOVE_FIXATIONS_FROM_GAMETES_HPP

#include <algorithm>
#include <cstdint>
#include <fwdpp/forward_types.hpp>

namespace fwdpp
{
    namespace ts
    {
        template <typename GenomeContainerType, typename MutationContainerType,
                  typename mutation_count_container>
        void
        remove_fixations_from_haploid_genomes(
            GenomeContainerType &haploid_genomes, const MutationContainerType &mutations,
            const mutation_count_container &mcounts,
            const mutation_count_container &mcounts_from_preserved_nodes,
            const fwdpp::uint_t twoN, const bool preserve_selected_fixations)
        /*! \brief Removed fixed variants from haploid_genomes
         *
         * This function should be called immediately after simplification and mutation counting.
         *
         * \param haploid_genomes A haploid_genome container
         * \param mutations A mutation container
         * \param mcounts A container stating how many times each element in \a mutations is present in the
         * currently-alive population
         * \param mcounts_from_preserved_nodes A container recording the counts of each mutation in ancient samples
         * \param twoN Twice the current population size
         * \param preserve_selected_fixations If true, do not mark selected fixations for recycling.
         *
         * A mutation is removed from a haploid_genome if one of the following conditions holds:
         * 1. The sum of \a mcounts and \a mcounts_from_preserved_nodes is zero.
         * 2. \a mcounts == \a twoN, mcounts_from_preserved_nodes is zero, 
         * and \a preserve_selected_fixations is false or the mutation is neutral.
         *
         * When \a preserve_selected_fixations is true, selected fixations are retained
         * in the population.  We do this because simulations of phenotypes (as opposed
         * to relative fitness) require tracking the contribution of fixation to 
         * trait values.
         *
         * \note When simulating a trait and never simulating neutral mutations, the most 
         * efficient thing is to skip calling this function entirely.
         *
         * \warning The value passed to \a preserve_selected_fixations needs to be 
         * coordinated with the \a preserve argument of flag_mutations_for_recycling.
         *
         * \todo Improve treatment of fixations by allowing for variants fixed in alive AND 
         * ancient samples to be detected.
         */
        {
            bool fixations_exist = false;
            for (std::size_t i = 0; !fixations_exist && i < mcounts.size(); ++i)
                {
                    if (mcounts[i] == twoN && mcounts_from_preserved_nodes[i] == 0)
                        {
                            fixations_exist = true;
                        }
                }
            if (fixations_exist)
                {
                    auto removal_criteria
                        = [&mutations, &mcounts, &mcounts_from_preserved_nodes, twoN,
                           preserve_selected_fixations](const fwdpp::uint_t key) {
                              return mcounts[key] == twoN
                                     && mcounts_from_preserved_nodes[key] == 0
                                     && (!preserve_selected_fixations
                                         || mutations[key].neutral);
                          };
                    for (auto &g : haploid_genomes)
                        {
                            if (g.n)
                                {
                                    auto itr = std::remove_if(g.mutations.begin(),
                                                              g.mutations.end(),
                                                              removal_criteria);
                                    g.mutations.erase(itr, g.mutations.end());
                                    itr = std::remove_if(g.smutations.begin(),
                                                         g.smutations.end(),
                                                         removal_criteria);
                                    g.smutations.erase(itr, g.smutations.end());
                                }
                        }
                }
        }
    } // namespace ts
} // namespace fwdpp

#endif
