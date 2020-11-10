#ifndef FWDPP_DATA_MATRIX_HPP_
#define FWDPP_DATA_MATRIX_HPP_

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
#include <fwdpp/poptypes/tags.hpp>
#include <gsl/gsl_matrix_char.h>

namespace fwdpp
{
    struct state_matrix
    /*! \brief Simplistic matrix representation of mutations
     * This type uses std::vector<std::int8_t> to hold a matrix
     * representing the genotypes for a set of diploids.
     *
     * For a haplotype matrix of n individuals, the data represent
     * 2n columns with a 0/1 encoding representing ancestral/derived.
     *
     * For a genotype matrix of n individuals, the data represent
     * n columns with a 0/1/2 encoding for the number of copies of the
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
     *
     * \ingroup samplingPops
     */
    {
        //! The state data
        std::vector<std::int8_t> data;
        //! Positions of variable sites
        std::vector<double> positions;
        //! Construct an empty object
        state_matrix() : data(), positions()
        {
        }
        //! Perfect-forwarding constructor
        template <typename T, typename P>
        state_matrix(T &&t, P &&p)
            : data(std::forward<T>(t)), positions(std::forward<P>(p))
        {
        }
    };

    struct data_matrix
    /*!
     * \brief Genotype or haplotype matrix.
     *
     * Data for neutral and selected variants,
     * respectively are stored as state_matrix objects.
     *
     * \ingroup samplingPops
     */
    {
        //! Data for neutral mutations.
        state_matrix neutral;
        //! Data for selected mutations.
        state_matrix selected;
        //! Locations of neutral mutations from mutation vector.  Same order
        //! as matrix row order
        std::vector<std::size_t> neutral_keys;
        //! Locations of selected mutations from mutation vector.  Same order
        //! as matrix row order
        std::vector<std::size_t> selected_keys;
        //! Number of columns in the matrix
        std::size_t ncol;
        explicit data_matrix(const std::size_t ncol_)
            /*!
             * Constructor
             *
             * \param ncol_ Number of columns in the matrix
             *
             * \version 0.7.0 
             * Changed from rows are sites to rows are individuals. Removed
             * default value of zero from constructor.
             */
            : neutral{}, selected{}, neutral_keys{}, selected_keys{}, ncol{ncol_}
        {
        }
    };
} // namespace fwdpp

// This header contains code re-used for
// implementing functions defined below.
#include "internal/data_matrix_details.hpp"

namespace fwdpp
{
    template <typename poptype>
    std::pair<std::vector<std::pair<std::size_t, uint_t>>,
              std::vector<std::pair<std::size_t, uint_t>>>
    mutation_keys(const poptype &pop, const std::vector<std::size_t> &individuals,
                  const bool include_neutral, const bool include_selected)
    /*!
     * For a sample defined by a set of diploids, obtain the keys corresponding
     * to all mutations in that sample.
     *
     * \param pop The population object
     * \param individuals The indexes of individuals in the sample
     * \param include_neutral If true, obtain keys for neutral sites
     * \param include_selected If true, obtain keys for selected sites
     *
     * \return A pair of vectors representing neutral and selected keys,
     * respectively.  The value
     * type of each vector is std::pair<std::size_t,fwdpp::uint_t>.  The first
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
     * the frequencies are summed, all using standard C++.
     *
     * \ingroup samplingPops
     */
    {
        return data_matrix_details::mutation_keys(
            pop.diploids, individuals, pop.haploid_genomes, pop.mcounts, include_neutral,
            include_selected, typename poptype::popmodel_t());
    }

    template <typename MutationContainerType>
    inline void
    sort_keys(const MutationContainerType &mutations,
              std::vector<std::pair<std::size_t, uint_t>> &keys)
    /*! \brief Sort keys by position
     *
     * Takes the data returned from mutation_keys and sorts
     * it by increasing mutation positions.
     *
     * \param mutations A mutation container
     * \param keys Returned from mutation_keys
     *
     * \ingroup samplingPops
     */
    {
        const auto comp = [&mutations](const std::pair<std::size_t, uint_t> &a,
                                       const std::pair<std::size_t, uint_t> &b) {
            return mutations[a.first].pos < mutations[b.first].pos;
        };
        std::sort(keys.begin(), keys.end(), comp);
    }

    template <typename F>
    inline void
    filter_keys(std::vector<std::pair<std::size_t, uint_t>> &keys, F f)
    /*! \brief Apply a filter to the keys
     *
     * Takes the data returned from mutation_keys and applies a 
     * filtering function.  The implementation is simply a wrapper
     * around the erase/remove idiom.
     *
     * \param keys 
     * \param f A function taking the value_type in \a keys and returning true if item should be removed
     */
    {
        keys.erase(std::remove_if(keys.begin(), keys.end(), f), keys.end());
    }

    template <typename poptype>
    data_matrix
    genotype_matrix(const poptype &pop, const std::vector<std::size_t> &individuals,
                    const std::vector<std::pair<std::size_t, uint_t>> &neutral_keys,
                    const std::vector<std::pair<std::size_t, uint_t>> &selected_keys)
    /*!
     * Calculate a fwdpp::data_matrix representing genotypes encoded as
     * 0,1, or 2 copies of the derived mutation.
     *
     * \param pop The population
     * \param individuals The indexes of individuals in \a pop forming the
     * sample.
     * \param neutral_keys See documentation of fwdpp::mutation_keys
     * \param selected_keys See documentation of fwdpp::mutation_keys
     *
     * \return fwdpp::data_matrix
     *
     * \ingroup samplingPops
     */
    {
        return data_matrix_details::fill_matrix(
            pop, individuals, neutral_keys, selected_keys,
            data_matrix_details::matrix_type::genotype);
    }

    template <typename poptype>
    data_matrix
    haplotype_matrix(const poptype &pop, const std::vector<std::size_t> &individuals,
                     const std::vector<std::pair<std::size_t, uint_t>> &neutral_keys,
                     const std::vector<std::pair<std::size_t, uint_t>> &selected_keys)
    /*!
     * Calculate a fwdpp::data_matrix representing haplotypes encoded as
     * 0 or 1 copies of the derived mutation.
     *
     * \param pop The population
     * \param individuals The indexes of individuals in \a pop forming the
     * sample.
     * \param neutral_keys See documentation of fwdpp::mutation_keys
     * \param selected_keys See documentation of fwdpp::mutation_keys
     *
     * \return fwdpp::data_matrix
     *
     * \ingroup samplingPops
     */
    {
        return data_matrix_details::fill_matrix(
            pop, individuals, neutral_keys, selected_keys,
            data_matrix_details::matrix_type::haplotype);
    }

    inline std::pair<std::vector<std::uint32_t>, std::vector<std::uint32_t>>
    row_sums(const data_matrix &m)
    /*!
     * Calculate the row sums of a fwdpp::data_matrix
     *
     * \return A pair of vectors of unsigned integers representing row sums
     * for neutral and selected sites in the matrix, respectively.
     *
     * \ingroup samplingPops
     */
    {
        return std::make_pair(
            data_matrix_details::row_col_sums_details(
                m.neutral.data, m.neutral.positions.size(), m.ncol, true),
            data_matrix_details::row_col_sums_details(
                m.selected.data, m.selected.positions.size(), m.ncol, true));
    }

    inline std::pair<std::vector<std::uint32_t>, std::vector<std::uint32_t>>
    col_sums(const data_matrix &m)
    /*!
     * Calculate the column sums of a fwdpp::data_matrix
     *
     * \return A pair of vectors of unsigned integers representing column
     * sums
     * for neutral and selected sites in the matrix, respectively.
     *
     * \ingroup samplingPops
     */
    {
        return std::make_pair(
            data_matrix_details::row_col_sums_details(
                m.neutral.data, m.neutral.positions.size(), m.ncol, false),
            data_matrix_details::row_col_sums_details(
                m.selected.data, m.selected.positions.size(), m.ncol, false));
    }
} // namespace fwdpp

#endif
