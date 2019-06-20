#include <fwdpp/poptypes/tags.hpp>

namespace fwdpp
{
    namespace fwdpp_internal
    {
        template <typename poptype>
        void
        zero_out_haploid_genomes(poptype &pop, fwdpp::poptypes::DIPLOID_TAG)
        {
            for (auto &dip : pop.diploids)
                {
                    pop.haploid_genomes[dip.first].n
                        = pop.haploid_genomes[dip.second].n = 0;
                }
        }
    } // namespace fwdpp_internal
} // namespace fwdpp
