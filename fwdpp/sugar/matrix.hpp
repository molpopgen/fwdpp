#ifndef FWDPP_MATRIX_HPP_
#define FWDPP_MATRIX_HPP_

#include <utility>
#include <vector>
#include <type_traits>
#include <cstdint>
#include <algorithm>
#include <stdexcept>
#include <unordered_set>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/type_traits.hpp>
#include <fwdpp/sugar/poptypes/tags.hpp>

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
        data_matrix(const std::size_t nrow_)
            : neutral{}, selected{}, neutral_positions{}, selected_positions{},
              neutral_popfreq{}, selected_popfreq{}, nrow{ nrow_ }
        {
        }
    };

    namespace data_matrix_details
    {
        enum class matrix_type
        //! Flag for which type of data matrix to fill
        {
            genotype,
            haplotype
        };

        struct matrix_helper
        //! Holds the data needed for generating data matrix
        {
            std::vector<std::size_t> neutral_keys, selected_keys;
            std::vector<std::int8_t> neutral_row, neutral_row2, selected_row,
                selected_row2;
            matrix_helper(const std::vector<std::size_t> &nk,
                          const std::vector<std::size_t> &sk)
                : neutral_keys{ nk }, selected_keys{ sk },
                  neutral_row(std::vector<std::int8_t>(nk.size(), 0)),
                  neutral_row2(std::vector<std::int8_t>(nk.size(), 0)),
                  selected_row(std::vector<std::int8_t>(sk.size(), 0)),
                  selected_row2(std::vector<std::int8_t>(sk.size(), 0))
            {
            }
            void
            zero()
            //! Refill temporary rows with zeros
            {
                std::fill(neutral_row.begin(), neutral_row.end(), 0);
                std::fill(neutral_row2.begin(), neutral_row2.end(), 0);
                std::fill(selected_row.begin(), selected_row.end(), 0);
                std::fill(selected_row2.begin(), selected_row2.end(), 0);
            }
        };

        template <typename gamete_t>
        void
        update_mutation_keys(std::unordered_set<std::size_t> &keys,
                             const typename gamete_t::mutation_container &a,
                             const std::vector<uint_t> &mcounts)
        {
            for (auto &&ai : a)
                {
                    if (mcounts[ai])
                        {
                            keys.insert(ai);
                        }
                }
        }

        template <typename dipvector_t, typename gcont_t>
        std::pair<std::vector<std::size_t>, std::vector<std::size_t>>
        mutation_keys(const dipvector_t &diploids,
                      const std::vector<std::size_t> &individuals,
                      const gcont_t &gametes,
                      const std::vector<uint_t> &mcounts,
                      const bool include_neutral, const bool include_selected,
                      const std::size_t, sugar::SINGLEPOP_TAG)
        {
            std::unordered_set<std::size_t> n, s;
            using gamete_t = typename gcont_t::value_type;
            for (auto &&ind : individuals)
                {
                    auto &dip = diploids[ind];
                    if (include_neutral)
                        {
                            update_mutation_keys<gamete_t>(
                                n, gametes[dip.first].mutations, mcounts);
                            update_mutation_keys<gamete_t>(
                                n, gametes[dip.second].mutations, mcounts);
                        }
                    if (include_selected)
                        {
                            update_mutation_keys<gamete_t>(
                                s, gametes[dip.first].smutations, mcounts);
                            update_mutation_keys<gamete_t>(
                                s, gametes[dip.second].smutations, mcounts);
                        }
                }
            return std::make_pair(
                std::vector<std::size_t>(n.begin(), n.end()),
                std::vector<std::size_t>(s.begin(), s.end()));
        }

        template <typename dipvector_t, typename gcont_t>
        std::pair<std::vector<std::size_t>, std::vector<std::size_t>>
        mutation_keys(const dipvector_t &diploids,
                      const std::vector<std::size_t> &individuals,
                      const gcont_t &gametes,
                      const std::vector<uint_t> &mcounts,
                      const bool include_neutral, const bool include_selected,
                      const std::size_t deme, sugar::METAPOP_TAG)
        {
            if (deme >= diploids.size())
                {
                    throw std::out_of_range("deme index out of range, "
                                            + std::string(__FILE__) + ' '
                                            + std::to_string(__LINE__));
                }
            // The i-th deme is a "singlepop", so we can use that function
            return mutation_keys(diploids[deme], individuals, gametes, mcounts,
                                 include_neutral, include_selected, deme,
                                 sugar::SINGLEPOP_TAG());
        }

        template <typename dipvector_t, typename gcont_t>
        std::pair<std::vector<std::size_t>, std::vector<std::size_t>>
        mutation_keys(const dipvector_t &diploids,
                      const std::vector<std::size_t> &individuals,
                      const gcont_t &gametes,
                      const std::vector<uint_t> &mcounts,
                      const bool include_neutral, const bool include_selected,
                      const std::size_t, sugar::MULTILOCPOP_TAG)
        {
            std::unordered_set<std::size_t> n, s;
            using gamete_t = typename gcont_t::value_type;
            for (auto &&ind : individuals)
                {
                    auto &dip = diploids[ind];
                    for (auto &&locus : dip)
                        {
                            if (include_neutral)
                                {
                                    update_mutation_keys<gamete_t>(
                                        n, gametes[locus.first].mutations,
                                        mcounts);
                                    update_mutation_keys<gamete_t>(
                                        n, gametes[locus.second].mutations,
                                        mcounts);
                                }
                            if (include_selected)
                                {
                                    update_mutation_keys<gamete_t>(
                                        s, gametes[locus.first].smutations,
                                        mcounts);
                                    update_mutation_keys<gamete_t>(
                                        s, gametes[locus.second].smutations,
                                        mcounts);
                                }
                        }
                }
            return std::make_pair(
                std::vector<std::size_t>(n.begin(), n.end()),
                std::vector<std::size_t>(s.begin(), s.end()));
        }

        inline void
        update_row(std::vector<std::int8_t> &v,
                   const std::vector<KTfwd::uint_t> &mut_keys,
                   const std::vector<std::size_t> &indexes)
        {
            if (v.size() != indexes.size())
                {
                    throw std::runtime_error("vector sizes do not match");
                }
            for (auto &&mk : mut_keys)
                {
                    auto i = std::find(indexes.begin(), indexes.end(),
                                       static_cast<std::size_t>(mk));
                    if (i != indexes.end()) // i may equal indexes.end iff mk
                                            // refers to a fixation
                        {
                            auto idx = std::distance(indexes.begin(), i);
                            if (static_cast<std::size_t>(idx) >= v.size())
                                {
                                    throw std::runtime_error(
                                        "idx >= v.size()");
                                }
                            v[static_cast<std::size_t>(idx)] += 1.0;
                        }
                }
        }

        template <typename gamete_t>
        inline void
        update_row_common(const gamete_t &g1, const gamete_t &g2,
                          matrix_helper &h)
        {
            update_row(h.neutral_row, g1.mutations, h.neutral_keys);
            update_row(h.neutral_row2, g2.mutations, h.neutral_keys);
            update_row(h.selected_row, g1.smutations, h.selected_keys);
            update_row(h.selected_row2, g2.smutations, h.selected_keys);
        }

        inline void
        fill_matrix_with_rows(data_matrix &m, matrix_helper &h,
                              const matrix_type mtype)
        {
            if (mtype == matrix_type::haplotype)
                {
                    m.neutral.insert(m.neutral.end(), h.neutral_row.begin(),
                                     h.neutral_row.end());
                    m.neutral.insert(m.neutral.end(), h.neutral_row2.begin(),
                                     h.neutral_row2.end());
                    m.selected.insert(m.selected.end(), h.selected_row.begin(),
                                      h.selected_row.end());
                    m.selected.insert(m.selected.end(),
                                      h.selected_row2.begin(),
                                      h.selected_row2.end());
                }
            else if (mtype == matrix_type::genotype)
                {
                    std::transform(h.neutral_row2.begin(),
                                   h.neutral_row2.end(), h.neutral_row.begin(),
                                   h.neutral_row.begin(),
                                   std::plus<std::int8_t>());
                    std::transform(
                        h.selected_row2.begin(), h.selected_row2.end(),
                        h.selected_row.begin(), h.selected_row.begin(),
                        std::plus<std::int8_t>());
                    m.neutral.insert(m.neutral.end(), h.neutral_row.begin(),
                                     h.neutral_row.end());
                    m.selected.insert(m.selected.end(), h.selected_row.begin(),
                                      h.selected_row.end());
                }
            h.zero();
        }

        template <typename poptype>
        void
        fill_matrix(const poptype &pop, data_matrix &m,
                    const std::vector<std::size_t> &individuals,
                    const std::vector<std::size_t> &neutral_keys,
                    const std::vector<std::size_t> &selected_keys,
                    const std::size_t, sugar::SINGLEPOP_TAG, matrix_type mtype)
        {
            matrix_helper h(neutral_keys, selected_keys);
            for (auto &&ind : individuals)
                {
                    if (ind >= pop.diploids.size())
                        {
                            throw std::out_of_range(
                                "individual index out of range");
                        }
                    auto &dip = pop.diploids[ind];
                    update_row_common(pop.gametes[dip.first],
                                      pop.gametes[dip.second], h);
                    fill_matrix_with_rows(m, h, mtype);
                }
            // fill out other data fields
            for (auto &&i : neutral_keys)
                {
                    m.neutral_positions.push_back(pop.mutations[i].pos);
                    m.neutral_popfreq.push_back(
                        static_cast<double>(pop.mcounts[i])
                        / static_cast<double>(2 * pop.diploids.size()));
                }
            for (auto &&i : selected_keys)
                {
                    m.selected_positions.push_back(pop.mutations[i].pos);
                    m.selected_popfreq.push_back(
                        static_cast<double>(pop.mcounts[i])
                        / static_cast<double>(2 * pop.diploids.size()));
                }
        }

        template <typename poptype>
        void
        fill_matrix(const poptype &pop, data_matrix &m,
                    const std::vector<std::size_t> &individuals,
                    const std::vector<std::size_t> &neutral_keys,
                    const std::vector<std::size_t> &selected_keys,
                    const std::size_t deme, sugar::METAPOP_TAG,
                    matrix_type mtype)
        {
            matrix_helper h(neutral_keys, selected_keys);
            for (auto &&ind : individuals)
                {
                    if (ind >= pop.diploids[deme].size())
                        {
                            throw std::out_of_range(
                                "individual index out of range");
                        }
                    auto &dip = pop.diploids[deme][ind];
                    update_row_common(pop.gametes[dip.first],
                                      pop.gametes[dip.second], h);
                    fill_matrix_with_rows(m, h, mtype);
                }
            // fill out other data fields
            for (auto &&i : neutral_keys)
                {
                    m.neutral_positions.push_back(pop.mutations[i].pos);
                    m.neutral_popfreq.push_back(
                        static_cast<double>(pop.mcounts[i])
                        / static_cast<double>(2 * pop.diploids[deme].size()));
                }
            for (auto &&i : selected_keys)
                {
                    m.selected_positions.push_back(pop.mutations[i].pos);
                    m.selected_popfreq.push_back(
                        static_cast<double>(pop.mcounts[i])
                        / static_cast<double>(2 * pop.diploids[deme].size()));
                }
        }
        template <typename poptype>
        void
        fill_matrix(const poptype &pop, data_matrix &m,
                    const std::vector<std::size_t> &individuals,
                    const std::vector<std::size_t> &neutral_keys,
                    const std::vector<std::size_t> &selected_keys,
                    const std::size_t, sugar::MULTILOCPOP_TAG,
                    matrix_type mtype)
        {
            matrix_helper h(neutral_keys, selected_keys);
            for (auto &&ind : individuals)
                {
                    if (ind >= pop.diploids.size())
                        {
                            throw std::out_of_range(
                                "individual index out of range");
                        }
                    auto &dip = pop.diploids[ind];
                    for (auto &&locus : dip)
                        {
                            update_row_common(pop.gametes[locus.first],
                                              pop.gametes[dip.second], h);
                        }
                    fill_matrix_with_rows(m, h, mtype);
                }
            // fill out other data fields
            for (auto &&i : neutral_keys)
                {
                    m.neutral_positions.push_back(pop.mutations[i].pos);
                    m.neutral_popfreq.push_back(
                        static_cast<double>(pop.mcounts[i])
                        / static_cast<double>(2 * pop.diploids.size()));
                }
            for (auto &&i : selected_keys)
                {
                    m.selected_positions.push_back(pop.mutations[i].pos);
                    m.selected_popfreq.push_back(
                        static_cast<double>(pop.mcounts[i])
                        / static_cast<double>(2 * pop.diploids.size()));
                }
        }

        template <typename poptype>
        data_matrix
        fill_matrix(const poptype &pop,
                    const std::vector<std::size_t> &individuals,
                    const std::vector<std::size_t> &neutral_keys,
                    const std::vector<std::size_t> &selected_keys,
                    const std::size_t deme, const matrix_type mtype)
        {
            data_matrix rv((mtype == matrix_type::genotype)
                               ? individuals.size()
                               : 2 * individuals.size());
            // dispatch details out depending on population type
            fill_matrix(pop, rv, individuals, neutral_keys, selected_keys,
                        deme, typename poptype::popmodel_t(), mtype);
            return rv;
        }
    }

    template <typename poptype>
    std::pair<std::vector<std::size_t>, std::vector<std::size_t>>
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
    genotype_matrix(const poptype &pop,
                    const std::vector<std::size_t> &individuals,
                    const std::vector<std::size_t> &neutral_keys,
                    const std::vector<std::size_t> &selected_keys,
                    const std::size_t deme = 0)
    {
        return data_matrix_details::fill_matrix(
            pop, individuals, neutral_keys, selected_keys, deme,
            data_matrix_details::matrix_type::genotype);
    }

    template <typename poptype>
    data_matrix
    haplotype_matrix(const poptype &pop,
                     const std::vector<std::size_t> &individuals,
                     const std::vector<std::size_t> &neutral_keys,
                     const std::vector<std::size_t> &selected_keys,
                     const std::size_t deme = 0)
    {
        return data_matrix_details::fill_matrix(
            pop, individuals, neutral_keys, selected_keys, deme,
            data_matrix_details::matrix_type::haplotype);
    }
}

#endif
