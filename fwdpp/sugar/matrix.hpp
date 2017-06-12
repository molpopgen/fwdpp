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
#include <gsl/gsl_matrix_char.h>

namespace KTfwd
{
    struct data_matrix
    /*!
     * \brief Genotype or haplotype matrix.
     *
     * This type uses std::vector<std::int8_t> to hold a matrix
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
        std::vector<std::int8_t> neutral;
        //! Data for selected mutations.
        std::vector<std::int8_t> selected;
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
        data_matrix(const std::size_t nrow_ = 0)
            /*!
             * Constructor
             *
             * \param nrow_ Number of rows in the matrix
             *
             * The default value of nrow_ = 0 is so that instances of this
             * type can be easily stack-allocated.  Intended use cases are
             * Cython,
             * Rcpp, or other systems for "wrapping" C++.  This type will
             * rely on compiler-generated copy/move constructors, and therefore
             * such
             * systems should rely on copy elision to make sure that
             * data_matrix::nrow
             * is set correctly.
             */
            : neutral{},
              selected{},
              neutral_positions{},
              selected_positions{},
              neutral_popfreq{},
              selected_popfreq{},
              nrow{ nrow_ }
        {
        }
    };
}

// This header contains code re-used for
// implementing functions defined below.
#include "matrix_details.hpp"

namespace KTfwd
{
    template <typename poptype>
    std::pair<std::vector<std::pair<std::size_t, uint_t>>,
              std::vector<std::pair<std::size_t, uint_t>>>
    mutation_keys(const poptype &pop,
                  const std::vector<std::size_t> &individuals,
                  const bool include_neutral, const bool include_selected,
                  const std::size_t deme = 0)
    /*!
     * For a sample defined by a set of diploids, obtain the keys corresponding
     * to all mutations in that sample.
     *
     * \param pop The population object
     * \param individuals The indexes of individuals in the sample
     * \param include_neutral If true, obtain keys for neutral sites
     * \param include_selected If true, obtain keys for selected sites
     * \param deme If poptype represents a metapopulation object, the
     * individuals come from this deme.
     *
     * \return A pair of vectors representing neutral and selected keys,
     * respectively.  The value
     * type of each vector is std::pair<std::size_t,KTfwd::uint_t>.  The first
     * value is the mutation
     * key, and the second is its frequency in the sample.
     *
     * Several comments are required:
     *
     * 1. Return values are unsorted with respect to anything meaningful.  If
     * you wish to sort on position,
     * etc., do so yourself.
     * 2. In general, any manipulation of the keys is possible.  For example,
     * removing mutations that are fixed
     * in the sample, or whose frequency is or is not within some desired
     * range, is easily doable via the
     * "erase/remove" idiom.
     * 3. Keys from multiple samples can be merged to form new vectors where
     * the key elements are unique and
     * the frequencies are summed, all using standard C++.  Example use case is
     * combining samples from
     * multiple demes.
     */
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
    /*!
     * Calculate a KTfwd::data_matrix representing genotypes encoded as
     * 0,1, or 2 copies of the derived mutation.
     *
     * \param pop The population
     * \param individuals The indexes of individuals in \a pop forming the
     * sample.
     * \param neutral_keys See documentation of KTfwd::mutation_keys
     * \param selected_keys See documentation of KTfwd::mutation_keys
     * \param deme If pop is a metapopulation, this is the deme containing the
     * sample
     *
     * \return KTfwd::data_matrix
     *
     * \note Return values representing samples from different demes must be
     * combined
     * by the user.
     */
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
    /*!
     * Calculate a KTfwd::data_matrix representing haplotypes encoded as
     * 0 or 1 copies of the derived mutation.
     *
     * \param pop The population
     * \param individuals The indexes of individuals in \a pop forming the
     * sample.
     * \param neutral_keys See documentation of KTfwd::mutation_keys
     * \param selected_keys See documentation of KTfwd::mutation_keys
     * \param deme If pop is a metapopulation, this is the deme containing the
     * sample
     *
     * \return KTfwd::data_matrix
     *
     * \note Return values representing samples from different demes must be
     * combined
     * by the user.
     */
    {
        return data_matrix_details::fill_matrix(
            pop, individuals, neutral_keys, selected_keys, deme,
            data_matrix_details::matrix_type::haplotype);
    }

    inline std::pair<std::vector<std::uint32_t>, std::vector<std::uint32_t>>
    row_sums(const data_matrix &m)
    /*!
     * Calculate the row sums of a KTfwd::data_matrix
     *
     * \return A pair of vectors of unsigned integers representing row sums
     * for neutral and selected sites in the matrix, respectively.
     */
    {
        return std::make_pair(
            data_matrix_details::row_col_sums_details(
                m.neutral, m.nrow, m.neutral_positions.size(), true),
            data_matrix_details::row_col_sums_details(
                m.selected, m.nrow, m.selected_positions.size(), true));
    }

    inline std::pair<std::vector<std::uint32_t>, std::vector<std::uint32_t>>
    col_sums(const data_matrix &m)
    /*!
         * Calculate the column sums of a KTfwd::data_matrix
         *
         * \return A pair of vectors of unsigned integers representing column
     * sums
         * for neutral and selected sites in the matrix, respectively.
         */
    {
        return std::make_pair(
            data_matrix_details::row_col_sums_details(
                m.neutral, m.nrow, m.neutral_positions.size(), false),
            data_matrix_details::row_col_sums_details(
                m.selected, m.nrow, m.selected_positions.size(), false));
    }
}

#endif
