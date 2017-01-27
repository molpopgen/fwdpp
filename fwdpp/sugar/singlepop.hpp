#ifndef __FWDPP_SUGAR_SINGLEPOP_HPP__
#define __FWDPP_SUGAR_SINGLEPOP_HPP__

#include <utility>
#include <vector>
#include <unordered_set>
#include <fwdpp/fwd_functional.hpp>
#include <fwdpp/sugar/poptypes/singlepop.hpp>

namespace KTfwd
{
    /*!
      \brief Single locus, single population object
      \ingroup sugar
    */
    template <typename mtype,
              typename diploid_t = std::pair<std::size_t, std::size_t>>
    using singlepop
        = sugar::singlepop<mtype, std::vector<mtype>, std::vector<gamete>,
                           std::vector<diploid_t>, std::vector<mtype>,
                           std::vector<uint_t>,
                           std::unordered_set<double, std::hash<double>,
                                              KTfwd::equal_eps>>;
}
#endif
