#ifndef __FWDPP_SUGAR_MULTILOC_HPP__
#define __FWDPP_SUGAR_MULTILOC_HPP__

#include <vector>
#include <unordered_set>
#include <fwdpp/sugar/poptypes/multiloc.hpp>
#include <fwdpp/fwd_functional.hpp>

namespace KTfwd
{
    /*!
      \brief Single population, multilocus simulation.
      See @ref md_md_sugar for rationale, etc.
      \ingroup sugar
    */
    template <typename mtype,
              typename diploid_t = std::pair<std::size_t, std::size_t>>
    using multiloc
        = sugar::multiloc<mtype, std::vector<mtype>, std::vector<gamete>,
                          std::vector<std::vector<diploid_t>>,
                          std::vector<mtype>, std::vector<uint_t>,
                          std::unordered_set<double, std::hash<double>,
                                             KTfwd::equal_eps>>;
}
#endif
