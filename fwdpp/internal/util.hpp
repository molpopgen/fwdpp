#include <fwdpp/poptypes/tags.hpp>

namespace fwdpp
{
    namespace fwdpp_internal
    {
        template <typename poptype>
        void
        zero_out_gametes(poptype &pop, fwdpp::poptypes::DIPLOID_TAG)
        {
            for (auto &dip : pop.diploids)
                {
                    pop.gametes[dip.first].n = pop.gametes[dip.second].n = 0;
                }
        }
    } // namespace fwdpp_internal
} // namespace fwdpp
