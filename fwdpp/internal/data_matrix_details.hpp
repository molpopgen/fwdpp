#ifndef FWDPP_DATA_MATRIX_DETAILS_HPP
#define FWDPP_DATA_MATRIX_DETAILS_HPP

#include <stdexcept>
#include <iterator>
#include <numeric>
#include <fwdpp/debug.hpp>

/*
 * This header is not meant to be included directly.
 */
namespace fwdpp
{
    namespace data_matrix_details
    {
        enum class matrix_type
        //! Flag for which type of data matrix to fill
        {
            genotype,
            haplotype
        };

        template <typename mutation_key_container>
        void
        update_mutation_keys(std::unordered_map<std::size_t, uint_t> &keys,
                             const mutation_key_container &a,
                             const std::vector<uint_t> &mcounts)
        {
            for (auto &&ai : a)
                {
                    if (mcounts[ai])
                        {
                            auto k = keys.find(ai);
                            if (k == keys.end())
                                {
                                    keys.emplace(ai, 1);
                                }
                            else
                                {
                                    k->second++;
                                }
                        }
                }
        }

        template <typename dipvector_t, typename gcont_t>
        std::pair<std::vector<std::pair<std::size_t, uint_t>>,
                  std::vector<std::pair<std::size_t, uint_t>>>
        mutation_keys(const dipvector_t &diploids,
                      const std::vector<std::size_t> &individuals,
                      const gcont_t &haploid_genomes,
                      const std::vector<uint_t> &mcounts,
                      const bool include_neutral, const bool include_selected,
                      poptypes::DIPLOID_TAG)
        {
            std::unordered_map<std::size_t, uint_t> n, s;
            for (auto &&ind : individuals)
                {
                    auto &dip = diploids[ind];
                    if (include_neutral)
                        {
                            update_mutation_keys(
                                n, haploid_genomes[dip.first].mutations,
                                mcounts);
                            update_mutation_keys(
                                n, haploid_genomes[dip.second].mutations,
                                mcounts);
                        }
                    if (include_selected)
                        {
                            update_mutation_keys(
                                s, haploid_genomes[dip.first].smutations,
                                mcounts);
                            update_mutation_keys(
                                s, haploid_genomes[dip.second].smutations,
                                mcounts);
                        }
                }
            return std::make_pair(std::vector<std::pair<std::size_t, uint_t>>(
                                      n.begin(), n.end()),
                                  std::vector<std::pair<std::size_t, uint_t>>(
                                      s.begin(), s.end()));
        }

        template <typename mcont_t, typename key_container>
        inline void
        update_pos(const mcont_t &mutations, const key_container &keys,
                   state_matrix &sm)
        {
            for (auto &key : keys)
                {
                    sm.positions.push_back(mutations[key.first].pos);
                }
        }

        template <typename mutation_key_container>
        void
        update_site(const mutation_key_container &first,
                    const mutation_key_container &second,
                    std::vector<std::int8_t> &site,
                    const std::pair<std::size_t, uint_t> &mutation_record,
                    const matrix_type mtype)
        {
            int onfirst
                = (std::find(first.begin(), first.end(), mutation_record.first)
                   != first.end());
            int onsecond = (std::find(second.begin(), second.end(),
                                      mutation_record.first)
                            != second.end());
            if (mtype == matrix_type::genotype)
                {
                    site.push_back(onfirst + onsecond);
                }
            else
                {
                    site.push_back(onfirst);
                    site.push_back(onsecond);
                }
        }

        template <typename poptype>
        void
        fill_matrix(
            const poptype &pop, data_matrix &m,
            const std::vector<std::size_t> &individuals,
            const std::vector<std::pair<std::size_t, uint_t>> &neutral_keys,
            const std::vector<std::pair<std::size_t, uint_t>> &selected_keys,
            poptypes::DIPLOID_TAG, matrix_type mtype)
        {
            for (auto &&mkey : neutral_keys)
                {
                    for (auto &ind : individuals)
                        {
                            update_site(
                                pop.haploid_genomes[pop.diploids[ind].first]
                                    .mutations,
                                pop.haploid_genomes[pop.diploids[ind].second]
                                    .mutations,
                                m.neutral.data, mkey, mtype);
                        }
                    m.neutral_keys.push_back(mkey.first);
                }
            for (auto &&mkey : selected_keys)
                {
                    for (auto &ind : individuals)
                        {
                            update_site(
                                pop.haploid_genomes[pop.diploids[ind].first]
                                    .smutations,
                                pop.haploid_genomes[pop.diploids[ind].second]
                                    .smutations,
                                m.selected.data, mkey, mtype);
                        }
                    m.selected_keys.push_back(mkey.first);
                }
            // fill out other data fields
            update_pos(pop.mutations, neutral_keys, m.neutral);
            update_pos(pop.mutations, selected_keys, m.selected);
        }

        template <typename poptype>
        data_matrix
        fill_matrix(
            const poptype &pop, const std::vector<std::size_t> &individuals,
            const std::vector<std::pair<std::size_t, uint_t>> &neutral_keys,
            const std::vector<std::pair<std::size_t, uint_t>> &selected_keys,
            const matrix_type mtype)
        {
            data_matrix rv((mtype == matrix_type::genotype)
                               ? individuals.size()
                               : 2 * individuals.size());
            // dispatch details out depending on population type
            fill_matrix(pop, rv, individuals, neutral_keys, selected_keys,
                        typename poptype::popmodel_t(), mtype);
            return rv;
        }

        inline std::vector<std::uint32_t>
        row_col_sums_details(const std::vector<std::int8_t> &data,
                             const std::size_t nrow, const std::size_t ncol,
                             const bool is_row_sums)
        {
            std::vector<std::uint32_t> rv;
            if (!data.empty())
                {
                    auto v = gsl_matrix_char_const_view_array(
                        reinterpret_cast<const char *>(data.data()), nrow,
                        ncol);
                    const std::size_t X
                        = (is_row_sums) ? v.matrix.size1 : v.matrix.size2;
                    for (std::size_t rc = 0; rc < X; ++rc)
                        {
                            gsl_vector_char_const_view view
                                = (is_row_sums)
                                      ? gsl_matrix_char_const_row(&v.matrix,
                                                                  rc)
                                      : gsl_matrix_char_const_column(&v.matrix,
                                                                     rc);
                            unsigned sum = 0;
                            for (std::size_t i = 0; i < view.vector.size; ++i)
                                {
                                    sum += gsl_vector_char_get(&view.vector,
                                                               i);
                                }
                            rv.push_back(sum);
                        }
                }
            return rv;
        }
    } // namespace data_matrix_details
} // namespace fwdpp

#endif
