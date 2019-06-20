#ifndef __SAMPLING_FUNCTIONS_HPP__
#define __SAMPLING_FUNCTIONS_HPP__

#include <vector>
#include <map>
#include <functional>
#include <algorithm>
#include <utility>
#include <string>
#include <type_traits>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fwdpp/type_traits.hpp>
#include <fwdpp/data_matrix.hpp>
#include "internal/sampling_functions_details.hpp"

/*! @defgroup samplingPops Functions related to taking samples from simulated
  populations
*/

namespace fwdpp
{
    template <typename poptype>
    data_matrix
    sample_individuals(const poptype &pop,
                       const std::vector<std::size_t> &individuals,
                       const bool include_neutral, const bool include_selected,
                       const bool remove_fixed)
    /*!
     * \brief Create a fwdpp::data_matrix for a set of individuals.
     *
     * \param pop A population
     * \param individuals indexes of individuals in \a pop
     * \param include_neutral If true, populate fwdpp::data_matrix::neutral
     * \param include_selected If true, populate fwdpp::data_matrix::selected
     * \param remove_fixed If true, remove variants that are fixed in the sample.
     *
     * \return fwdpp::data_matrix
     *
     * \note The return value is a haplotype matrix.
     *
     * \ingroup samplingPops
     */
    {
        auto keys = fwdpp_internal::generate_filter_sort_keys(
            pop, individuals, include_neutral, include_selected, remove_fixed);
        return haplotype_matrix(pop, individuals, keys.first, keys.second);
    }
}
#endif
