#include <stdexcept>
#include <iterator>
#include <numeric>
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
                      const gcont_t &gametes,
                      const std::vector<uint_t> &mcounts,
                      const bool include_neutral, const bool include_selected,
                      sugar::SINGLELOC_TAG)
        {
            std::unordered_map<std::size_t, uint_t> n, s;
            for (auto &&ind : individuals)
                {
                    auto &dip = diploids[ind];
                    if (include_neutral)
                        {
                            update_mutation_keys(
                                n, gametes[dip.first].mutations, mcounts);
                            update_mutation_keys(
                                n, gametes[dip.second].mutations, mcounts);
                        }
                    if (include_selected)
                        {
                            update_mutation_keys(
                                s, gametes[dip.first].smutations, mcounts);
                            update_mutation_keys(
                                s, gametes[dip.second].smutations, mcounts);
                        }
                }
            return std::make_pair(std::vector<std::pair<std::size_t, uint_t>>(
                                      n.begin(), n.end()),
                                  std::vector<std::pair<std::size_t, uint_t>>(
                                      s.begin(), s.end()));
        }

        template <typename dipvector_t, typename gcont_t>
        std::pair<std::vector<std::pair<std::size_t, uint_t>>,
                  std::vector<std::pair<std::size_t, uint_t>>>
        mutation_keys(const dipvector_t &diploids,
                      const std::vector<std::size_t> &individuals,
                      const gcont_t &gametes,
                      const std::vector<uint_t> &mcounts,
                      const bool include_neutral, const bool include_selected,
                      sugar::MULTILOC_TAG)
        {
            std::unordered_map<std::size_t, uint_t> n, s;
            for (auto &&ind : individuals)
                {
                    auto &dip = diploids[ind];
                    for (auto &&locus : dip)
                        {
                            if (include_neutral)
                                {
                                    update_mutation_keys(
                                        n, gametes[locus.first].mutations,
                                        mcounts);
                                    update_mutation_keys(
                                        n, gametes[locus.second].mutations,
                                        mcounts);
                                }
                            if (include_selected)
                                {
                                    update_mutation_keys(
                                        s, gametes[locus.first].smutations,
                                        mcounts);
                                    update_mutation_keys(
                                        s, gametes[locus.second].smutations,
                                        mcounts);
                                }
                        }
                }
            return std::make_pair(std::vector<std::pair<std::size_t, uint_t>>(
                                      std::make_move_iterator(n.begin()),
                                      std::make_move_iterator(n.end())),
                                  std::vector<std::pair<std::size_t, uint_t>>(
                                      std::make_move_iterator(s.begin()),
                                      std::make_move_iterator(s.end())));
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
            sugar::SINGLELOC_TAG, matrix_type mtype)
        {
            for (auto &&mkey : neutral_keys)
                {
                    for (auto &ind : individuals)
                        {
                            update_site(
                                pop.gametes[pop.diploids[ind].first].mutations,
                                pop.gametes[pop.diploids[ind].second]
                                    .mutations,
                                m.neutral, mkey, mtype);
                        }
                }
            for (auto &&mkey : selected_keys)
                {
                    for (auto &ind : individuals)
                        {
                            update_site(pop.gametes[pop.diploids[ind].first]
                                            .smutations,
                                        pop.gametes[pop.diploids[ind].second]
                                            .smutations,
                                        m.selected, mkey, mtype);
                        }
                }
            // fill out other data fields
            for (auto &&i : neutral_keys)
                {
                    assert(pop.mutations[i.first].neutral);
                    m.neutral_positions.push_back(pop.mutations[i.first].pos);
                    m.neutral_popfreq.push_back(
                        static_cast<double>(pop.mcounts[i.first])
                        / static_cast<double>(2 * pop.diploids.size()));
                }
            for (auto &&i : selected_keys)
                {
                    assert(!pop.mutations[i.first].neutral);
                    m.selected_positions.push_back(pop.mutations[i.first].pos);
                    m.selected_popfreq.push_back(
                        static_cast<double>(pop.mcounts[i.first])
                        / static_cast<double>(2 * pop.diploids.size()));
                }
        }

        template <typename poptype>
        void
        fill_matrix(
            const poptype &pop, data_matrix &m,
            const std::vector<std::size_t> &individuals,
            const std::vector<std::pair<std::size_t, uint_t>> &neutral_keys,
            const std::vector<std::pair<std::size_t, uint_t>> &selected_keys,
            sugar::MULTILOC_TAG, matrix_type mtype)
        {
            for (auto &mkey : neutral_keys)
                {
                    assert(pop.mutations[mkey.first].neutral);
                    assert(pop.mcounts[mkey.first]);
                    //We need to find out what locus this mutation is in
                    double mpos = pop.mutations[mkey.first].pos;
                    std::size_t locus_index
                        = std::numeric_limits<std::size_t>::max();
                    for (std::size_t i = 0;
                         i < pop.locus_boundaries.size()
                         && locus_index
                                == std::numeric_limits<std::size_t>::max();
                         ++i)
                        {
                            if (mpos >= pop.locus_boundaries[i].first
                                && mpos < pop.locus_boundaries[i].second)
                                {
                                    locus_index = i;
                                }
                        }
                    if (locus_index == std::numeric_limits<std::size_t>::max())
                        {
                            throw std::runtime_error(
                                "mutation position not found in "
                                "pop.locus_boundaries");
                        }
                    auto current_size = m.neutral.size();
                    for (auto &ind : individuals)
                        {
                            auto locus = pop.diploids[ind][locus_index];
                            assert(pop.gametes[locus.first].n);
                            assert(pop.gametes[locus.second].n);
                            update_site(pop.gametes[locus.first].mutations,
                                        pop.gametes[locus.second].mutations,
                                        m.neutral, mkey, mtype);
                        }
                    if (std::accumulate(m.neutral.begin() + current_size,
                                        m.neutral.end(), 0)
                        == 0)
                        {
                            throw std::runtime_error(
                                "no variation found at site in this sample");
                        }
                }
            for (auto &mkey : selected_keys)
                {
                    assert(pop.mcounts[mkey.first]);
                    assert(!pop.mutations[mkey.first].neutral);
                    //We need to find out what locus this mutation is in
                    double mpos = pop.mutations[mkey.first].pos;
                    std::size_t locus_index
                        = std::numeric_limits<std::size_t>::max();
                    for (std::size_t i = 0;
                         i < pop.locus_boundaries.size()
                         && locus_index
                                == std::numeric_limits<std::size_t>::max();
                         ++i)
                        {
                            if (mpos >= pop.locus_boundaries[i].first
                                && mpos < pop.locus_boundaries[i].second)
                                {
                                    locus_index = i;
                                }
                        }
                    if (locus_index == std::numeric_limits<std::size_t>::max())
                        {
                            throw std::runtime_error(
                                "mutation position not found in "
                                "pop.locus_boundaries");
                        }
                    auto current_size = m.selected.size();
                    for (auto &ind : individuals)
                        {
                            auto locus = pop.diploids[ind][locus_index];
                            assert(pop.gametes[locus.first].n);
                            assert(pop.gametes[locus.second].n);
                            update_site(pop.gametes[locus.first].smutations,
                                        pop.gametes[locus.second].smutations,
                                        m.selected, mkey, mtype);
                        }
                    if (std::accumulate(m.selected.begin() + current_size,
                                        m.selected.end(), 0)
                        == 0)
                        {
                            throw std::runtime_error(
                                "no variation found at site in this sample");
                        }
                }
            // fill out other data fields
            for (auto &&i : neutral_keys)
                {
                    assert(pop.mutations[i.first].neutral);
                    m.neutral_positions.push_back(pop.mutations[i.first].pos);
                    m.neutral_popfreq.push_back(
                        static_cast<double>(pop.mcounts[i.first])
                        / static_cast<double>(2 * pop.diploids.size()));
                }
            for (auto &&i : selected_keys)
                {
                    assert(!pop.mutations[i.first].neutral);
                    m.selected_positions.push_back(pop.mutations[i.first].pos);
                    m.selected_popfreq.push_back(
                        static_cast<double>(pop.mcounts[i.first])
                        / static_cast<double>(2 * pop.diploids.size()));
                }
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
