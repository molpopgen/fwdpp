#ifndef __FWDPP_SUGAR_MULTILOC_HPP__
#define __FWDPP_SUGAR_MULTILOC_HPP__

#include <vector>
#include <unordered_map>
#include <fwdpp/poptypes/mlocuspop.hpp>
#include <fwdpp/fwd_functional.hpp>

namespace fwdpp
{
    /*!
      \brief Multilocus simulation.
    */
    template <typename mtype,
              typename diploid_t = std::pair<std::size_t, std::size_t>>
    using mlocuspop
        = poptypes::mlocuspop<mtype, std::vector<mtype>, std::vector<gamete>,
                           std::vector<std::vector<diploid_t>>,
                           std::vector<mtype>, std::vector<uint_t>,
						   // fwdpp 0.6.1 changed this from an unordered_set,
						   // in order to address a rare bug. See GitHub
						   // issue 130 for details.
                           std::unordered_multimap<double, uint_t>>;
}
#endif
