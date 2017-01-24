#ifndef FWDPP_INTERNAL_SAMPLE_DIPLOID_HELPERS
#define FWDPP_INTERNAL_SAMPLE_DIPLOID_HELPERS

#include <vector>

namespace KTfwd
{
    namespace fwdpp_internal
    {
        template <typename gcont_t, typename mcont_t>
        inline void
        process_gametes(const gcont_t &gametes, const mcont_t &mutations,
                        std::vector<uint_t> &mcounts)
        /*!
          For every non-extinct gamete, increment the count of its mutations
          using the frequency of the gamete.

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
            for (const auto &g : gametes)
                {
                    const auto n = g.n;
                    if (n) // only do this for extant gametes
                        {
                            for (const auto &m : g.mutations)
                                mcounts[m] += n;
                            for (const auto &m : g.smutations)
                                mcounts[m] += n;
                        }
                }
        }
    }
}

#endif
