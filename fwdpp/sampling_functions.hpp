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

    template <typename poptype>
    std::vector<data_matrix>
    sample_individuals_by_window(
        const poptype &pop, const std::vector<std::size_t> &individuals,
        const std::vector<std::pair<double, double>> window_boundaries,
        const bool include_neutral, const bool include_selected,
        const bool remove_fixed)
    /*!
     * \brief Create a vector fwdpp::data_matrix for a set of individuals.
     *
     * \param pop A population
     * \param individuals indexes of individuals in \a pop
     * \param window_boundaries [start,stop) positions of each window.
     * \param include_neutral If true, populate fwdpp::data_matrix::neutral
     * \param include_selected If true, populate fwdpp::data_matrix::selected
     * \param remove_fixed If true, remove variants that are fixed in the sample.
     *
     * \return vector of fwdpp::data_matrix
     *
     * \note Each fwdpp::data_matrix is in haplotype matrix layout.
     *
     * This function differs from sample_individuals in that a 
     * separate data_matrix is returned for each window.  The main 
     * use case envisioned is for sampling from multilocus populations,
     * where mlocuspop::locus_boundaries is passed in as \a window_boundaries.
     * Other applications are possible, however.
     *
     * \ingroup samplingPops
     */
    {
        auto keys = fwdpp_internal::generate_filter_sort_keys(
            pop, individuals, include_neutral, include_selected, remove_fixed);
        std::vector<data_matrix> rv;
        decltype(keys.first.begin()) nstart, nend, sstart, send;

        const auto lbf = [&pop](const std::pair<std::size_t, uint_t> &p,
                                const double pos) {
            return pop.mutations[p.first].pos < pos;
        };

        const auto ubf = [&pop](const double pos,
                                const std::pair<std::size_t, uint_t> &p) {
            return pos < pop.mutations[p.first].pos;
        };
        decltype(keys.first) nwk, swk;
        /* Note:
         * The binary searches work for the case
         * of an empty window because, for e.g.,
         * you end up with nstart == nend,
         * which has the effect of clearing
         * nwk's contents.
         *
         * This case is tested in the test suite.
         */
        for (const auto &b : window_boundaries)
            {
                nstart = std::lower_bound(keys.first.begin(), keys.first.end(),
                                          b.first, lbf);
                nend = std::upper_bound(nstart, keys.first.end(), b.second,
                                        ubf);
                sstart = std::lower_bound(keys.second.begin(),
                                          keys.second.end(), b.first, lbf);
                send = std::upper_bound(sstart, keys.second.end(), b.second,
                                        ubf);
                nwk.assign(nstart, nend);
                swk.assign(sstart, send);

                rv.emplace_back(haplotype_matrix(pop, individuals, nwk, swk));
            }
        return rv;
    }
}
#endif
