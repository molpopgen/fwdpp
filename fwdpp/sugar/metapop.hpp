#ifndef __FWDPP_SUGAR_METAPOP_HPP__
#define __FWDPP_SUGAR_METAPOP_HPP__

#include <fwdpp/sugar/poptypes/metapop.hpp>
#include <fwdpp/fwd_functional.hpp>
#include <vector>
#include <unordered_set>

namespace KTfwd
{
    /*!
      \brief Single locus metapopulation simulation object
      \ingroup sugar
    */
    template <typename mtype,
              typename diploid_t = std::pair<std::size_t, std::size_t>>
    using metapop
        = sugar::metapop<mtype, std::vector<mtype>, std::vector<gamete>,
                         std::vector<diploid_t>,
                         std::vector<std::vector<diploid_t>>,
                         std::vector<mtype>, std::vector<uint_t>,
                         std::unordered_set<double, std::hash<double>,
                                            KTfwd::equal_eps>>;
}
#endif
