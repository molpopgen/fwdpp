#ifndef FWDPP_INTERNAL_SAMPLE_DIPLOID_HELPERS
#define FWDPP_INTERNAL_SAMPLE_DIPLOID_HELPERS

#include <vector>
#include <fwdpp/fundamental_types/typedefs.hpp>

namespace fwdpp
{
    namespace fwdpp_internal
    {
        template <typename GenomeContainerType, typename MutationContainerType>
        inline void
        process_haploid_genomes(const GenomeContainerType &haploid_genomes,
                                const MutationContainerType &mutations,
                                std::vector<uint_t> &mcounts)
        /*!
          For every non-extinct haploid_genome, increment the count of its mutations
          using the frequency of the haploid_genome.

          This is usually the most expensive function call in a simulation.
        */
        {
            // zero out mcounts
            for (auto &mc : mcounts)
                mc = 0;
            if (mutations.size() > mcounts.size())
                {
                    mcounts.resize(mutations.size(), 0);
                }
            // update mutation counts
            for (const auto &g : haploid_genomes)
                {
                    const auto n = g.n;
                    if (n) // only do this for extant haploid_genomes
                        {
                            for (const auto &m : g.mutations)
                                mcounts[m] += n;
                            for (const auto &m : g.smutations)
                                mcounts[m] += n;
                        }
                }
        }
    } // namespace fwdpp_internal
} // namespace fwdpp

#endif
