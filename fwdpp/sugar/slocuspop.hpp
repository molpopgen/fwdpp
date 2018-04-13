#ifndef __FWDPP_SUGAR_SINGLEPOP_HPP__
#define __FWDPP_SUGAR_SINGLEPOP_HPP__

#include <utility>
#include <vector>
#include <unordered_set>
#include <fwdpp/fwd_functional.hpp>
#include <fwdpp/sugar/poptypes/slocuspop.hpp>

namespace fwdpp
{
    /*!
      \brief Single locus, single population object
      \example juvenile_migration.cc
      \example K_linked_regions_extensions.cc
      \example K_linked_regions_generalized_rec.cc
      \ingroup sugar
    */
    template <typename mtype,
              typename diploid_t = std::pair<std::size_t, std::size_t>>
    using slocuspop
        = sugar::slocuspop<mtype, std::vector<mtype>, std::vector<gamete>,
                           std::vector<diploid_t>, std::vector<mtype>,
                           std::vector<uint_t>,
                           std::unordered_set<double, std::hash<double>,
                                              fwdpp::equal_eps>>;
}
#endif
