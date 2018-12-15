#include <fwdpp/poptypes/tags.hpp>

namespace fwdpp
{
    namespace fwdpp_internal
    {
        template <typename poptype>
        void
        zero_out_gametes(poptype &pop, fwdpp::poptypes::SINGLELOC_TAG)
        {
            for (auto &dip : pop.diploids)
                {
                    pop.gametes[dip.first].n = pop.gametes[dip.second].n = 0;
                }
        }

        template <typename poptype>
        void
        zero_out_gametes(poptype &pop, fwdpp::poptypes::MULTILOC_TAG)
        {
            for (auto &dip : pop.diploids)
                {
                    for (auto &locus : dip)
                        {
                            pop.gametes[locus.first].n = 0;
                            pop.gametes[locus.second].n = 0;
                        }
                }
        }
    } // namespace fwdpp_internal
} // namespace fwdpp
