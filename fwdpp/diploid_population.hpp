#ifndef __FWDPP_SUGAR_DIPLOID_POPULATION_HPP__
#define __FWDPP_SUGAR_DIPLOID_POPULATION_HPP__

#include <utility>
#include <vector>
#include <unordered_map>
#include <fwdpp/fwd_functional.hpp>
#include <fwdpp/poptypes/diploid_population.hpp>

namespace fwdpp
{
    /*!
      \brief Single locus, single population object
      \example juvenile_migration.cc
      \example K_linked_regions_extensions.cc
      \example K_linked_regions_generalized_rec.cc
    */
    template <typename mtype,
              typename diploid_t = std::pair<std::size_t, std::size_t>>
    using diploid_population = poptypes::diploid_population<
        mtype, std::vector<mtype>, std::vector<haploid_genome>,
        std::vector<diploid_t>, std::vector<mtype>, std::vector<uint_t>,
        // fwdpp 0.6.1 changed this from an unordered_set,
        // in order to address a rare bug. See GitHub
        // issue 130 for details.
        std::unordered_multimap<double, std::uint32_t>>;
} // namespace fwdpp
#endif
