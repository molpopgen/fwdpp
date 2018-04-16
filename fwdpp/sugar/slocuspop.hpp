#ifndef __FWDPP_SUGAR_SINGLEPOP_HPP__
#define __FWDPP_SUGAR_SINGLEPOP_HPP__

#include <utility>
#include <vector>
#include <unordered_map>
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
						   // fwdpp 0.6.1 changed this from an unordered_set,
						   // in order to address a rare bug. See GitHub
						   // issue 130 for details.
                           std::unordered_multimap<double, std::uint32_t>>;
}
#endif
