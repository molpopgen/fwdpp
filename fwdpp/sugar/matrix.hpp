#ifndef FWDPP_MATRIX_HPP_
#define FWDPP_MATRIX_HPP_

#include <cassert>
#include <numeric>
#include <utility>
#include <vector>
#include <type_traits>
#include <cstdint>
#include <algorithm>
#include <stdexcept>
#include <set>
#include <unordered_map>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/type_traits.hpp>
#include <fwdpp/sugar/poptypes/tags.hpp>
#include <gsl/gsl_matrix_short.h>

namespace KTfwd
{
    struct data_matrix
    /*!
     * \brief Genotype or haplotype matrix.
     *
     * This type uses std::vector<short> to hold a matrix
     * representing the genotypes for a set of diploids.
     *
     * For a haplotype matrix of n individuals, the data represent
     * 2n rows with a 0/1 encoding representing ancetral/derived.
     *
     * For a genotype matrix of n individuals, the data represent
     * n rows with a 0/1/2 encoding for the number of copies of the
     * derived mutation.
     *
     * The data layout is row-major (aka "C-style") ordering,
     * facilitating compatibility with GSL matrix types,
     * NumPy 2D arrays, etc.  Note that GSL matrices may be
	 * constructed using gsl_matrix_view_array or
	 * gsl_matrix_const_view_array for cases where a matrix of
	 * only neutral or only selected mutations is needed.
	 *
	 * We use the 8-bit integer type to save space.  In practice,
	 * one may convert (via copy) to other types for operations like
	 * regression.
     *
     * \note This type is not constructed directly, but rather returned
     * by other functions.
     */
    {
        //! Data for neutral mutations.
        std::vector<short> neutral;
        //! Data for selected mutations.
        std::vector<short> selected;
        //! Positions of neutral mutations.  Same order as matrix column order
        std::vector<double> neutral_positions;
        //! Positions of selected mutations.  Same order as matrix column order
        std::vector<double> selected_positions;
        //! Frequencies of neutral mutations in entire population.  Same order
        //! as matrix column order
        std::vector<double> neutral_popfreq;
        //! Frequencies of selected mutations in entire population.  Same order
        //! as matrix column order
        std::vector<double> selected_popfreq;
        //! Number of rows in the matrix
        std::size_t nrow;
        data_matrix(const std::size_t nrow_)
            : neutral{}, selected{}, neutral_positions{}, selected_positions{},
              neutral_popfreq{}, selected_popfreq{}, nrow{ nrow_ }
        {
        }
    };
}

#include "matrix_details.hpp"

namespace KTfwd {
    template <typename poptype>
    std::pair<std::vector<std::pair<std::size_t, uint_t>>,
              std::vector<std::pair<std::size_t, uint_t>>>
    mutation_keys(const poptype &pop,
                  const std::vector<std::size_t> &individuals,
                  const bool include_neutral, const bool include_selected,
                  const std::size_t deme = 0)
    {
        return data_matrix_details::mutation_keys(
            pop.diploids, individuals, pop.gametes, pop.mcounts,
            include_neutral, include_selected, deme,
            typename poptype::popmodel_t());
    }

    template <typename poptype>
    data_matrix
    genotype_matrix(
        const poptype &pop, const std::vector<std::size_t> &individuals,
        const std::vector<std::pair<std::size_t, uint_t>> &neutral_keys,
        const std::vector<std::pair<std::size_t, uint_t>> &selected_keys,
        const std::size_t deme = 0)
    {
        return data_matrix_details::fill_matrix(
            pop, individuals, neutral_keys, selected_keys, deme,
            data_matrix_details::matrix_type::genotype);
    }

    template <typename poptype>
    data_matrix
    haplotype_matrix(
        const poptype &pop, const std::vector<std::size_t> &individuals,
        const std::vector<std::pair<std::size_t, uint_t>> &neutral_keys,
        const std::vector<std::pair<std::size_t, uint_t>> &selected_keys,
        const std::size_t deme = 0)
    {
        return data_matrix_details::fill_matrix(
            pop, individuals, neutral_keys, selected_keys, deme,
            data_matrix_details::matrix_type::haplotype);
    }

    inline std::pair<std::vector<std::uint32_t>, std::vector<std::uint32_t>>
    row_sums(const data_matrix &m)
    {
        return std::make_pair(
            data_matrix_details::row_col_sums_details(
                m.neutral, m.nrow, m.neutral_positions.size(), true),
            data_matrix_details::row_col_sums_details(
                m.selected, m.nrow, m.selected_positions.size(), true));
    }
    inline std::pair<std::vector<std::uint32_t>, std::vector<std::uint32_t>>
    col_sums(const data_matrix &m)
    {
        return std::make_pair(
            data_matrix_details::row_col_sums_details(
                m.neutral, m.nrow, m.neutral_positions.size(), false),
            data_matrix_details::row_col_sums_details(
                m.selected, m.nrow, m.selected_positions.size(), false));
    }
}

#endif
