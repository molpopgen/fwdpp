#include <iterator>
/*
 * This header is not meant to be included directly.
 */
namespace KTfwd
{
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
            std::vector<std::pair<std::size_t, uint_t>> neutral_keys,
                selected_keys;
            std::vector<std::int8_t> neutral_row, neutral_row2, selected_row,
                selected_row2;
            matrix_helper(
                const std::vector<std::pair<std::size_t, uint_t>> &nk,
                const std::vector<std::pair<std::size_t, uint_t>> &sk)
                : neutral_keys{ nk }, selected_keys{ sk },
                  neutral_row(std::vector<std::int8_t>(nk.size(), 0)),
                  neutral_row2(std::vector<std::int8_t>(nk.size(), 0)),
                  selected_row(std::vector<std::int8_t>(sk.size(), 0)),
                  selected_row2(std::vector<std::int8_t>(sk.size(), 0))
            {
            }
            void
            zero()
            //! Refill temporary rows with zerosconst std::vector<uint_t> &
            //! gamete_mut_keys,
            {
                std::fill(neutral_row.begin(), neutral_row.end(), 0);
                std::fill(neutral_row2.begin(), neutral_row2.end(), 0);
                std::fill(selected_row.begin(), selected_row.end(), 0);
                std::fill(selected_row2.begin(), selected_row2.end(), 0);
            }
        };

        template <typename gamete_t>
        void
        update_mutation_keys(std::unordered_map<std::size_t, uint_t> &keys,
                             const typename gamete_t::mutation_container &a,
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
                      const std::size_t, sugar::SINGLEPOP_TAG)
        {
            std::unordered_map<std::size_t, uint_t> n, s;
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
        std::pair<std::vector<std::pair<std::size_t, uint_t>>,
                  std::vector<std::pair<std::size_t, uint_t>>>
        mutation_keys(const dipvector_t &diploids,
                      const std::vector<std::size_t> &individuals,
                      const gcont_t &gametes,
                      const std::vector<uint_t> &mcounts,
                      const bool include_neutral, const bool include_selected,
                      const std::size_t, sugar::MULTILOCPOP_TAG)
        {
            std::unordered_map<std::size_t, uint_t> n, s;
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
            return std::make_pair(std::vector<std::pair<std::size_t, uint_t>>(
                                      std::make_move_iterator(n.begin()),
                                      std::make_move_iterator(n.end())),
                                  std::vector<std::pair<std::size_t, uint_t>>(
                                      std::make_move_iterator(s.begin()),
                                      std::make_move_iterator(s.end())));
        }

        inline void
        update_row(std::vector<std::int8_t> &v,
                   const std::vector<KTfwd::uint_t> &mut_keys,
                   const std::vector<std::pair<std::size_t, uint_t>> &indexes)
        {
            if (v.size() != indexes.size())
                {
                    throw std::runtime_error("vector sizes do not match");
                }
            for (auto &&mk : mut_keys)
                {
                    auto i = std::find_if(
                        indexes.begin(), indexes.end(),
                        [mk](const std::pair<std::size_t, uint_t> &p) {
                            return p.first == mk;
                        });
                    if (i != indexes.end()) // i may equal indexes.end iff mk
                        // refers to a fixation
                        {
                            auto idx = std::distance(indexes.begin(), i);
                            if (static_cast<std::size_t>(idx) >= v.size())
                                {
                                    throw std::runtime_error(
                                        "idx >= v.size()");
                                }
                            v[static_cast<std::size_t>(idx)]++;
                        }
                }
        }

#ifndef NDEBUG
        inline bool
        validate_rows(const std::vector<uint_t> &gamete_mut_keys,
                      const std::vector<std::pair<std::size_t, uint_t>> &keys,
                      const std::vector<std::int8_t> &row)
        //! check that row sums are ok.
        // We need this more expensive check in case keys are adjusted prior
        // to filling matrix.
        {
            std::set<std::size_t> gam(gamete_mut_keys.begin(),
                                      gamete_mut_keys.end());
            std::set<std::size_t> k;
            for (auto &&ki : keys)
                {
                    k.insert(ki.first);
                }
            std::vector<std::size_t> intersection;
            std::set_intersection(gam.begin(), gam.end(), k.begin(), k.end(),
                                  std::back_inserter(intersection));
            return intersection.size()
                   == std::accumulate(row.begin(), row.end(), 0.);
        }
#endif

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
                                   h.neutral_row.begin(), std::plus<std::int8_t>());
                    std::transform(h.selected_row2.begin(),
                                   h.selected_row2.end(),
                                   h.selected_row.begin(),
                                   h.selected_row.begin(), std::plus<std::int8_t>());
                    m.neutral.insert(m.neutral.end(), h.neutral_row.begin(),
                                     h.neutral_row.end());
                    m.selected.insert(m.selected.end(), h.selected_row.begin(),
                                      h.selected_row.end());
                }
            h.zero();
        }

        template <typename poptype>
        void
        fill_matrix(
            const poptype &pop, data_matrix &m,
            const std::vector<std::size_t> &individuals,
            const std::vector<std::pair<std::size_t, uint_t>> &neutral_keys,
            const std::vector<std::pair<std::size_t, uint_t>> &selected_keys,
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
                    assert(validate_rows(pop.gametes[dip.first].mutations,
                                         h.neutral_keys, h.neutral_row));
                    assert(validate_rows(pop.gametes[dip.second].mutations,
                                         h.neutral_keys, h.neutral_row2));
                    assert(validate_rows(pop.gametes[dip.first].smutations,
                                         h.selected_keys, h.selected_row));
                    assert(validate_rows(pop.gametes[dip.second].smutations,
                                         h.selected_keys, h.selected_row2));
                    fill_matrix_with_rows(m, h, mtype);
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
            const std::size_t deme, sugar::METAPOP_TAG, matrix_type mtype)
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
                    assert(validate_rows(pop.gametes[dip.first].mutations,
                                         h.neutral_keys, h.neutral_row));
                    assert(validate_rows(pop.gametes[dip.second].mutations,
                                         h.neutral_keys, h.neutral_row2));
                    assert(validate_rows(pop.gametes[dip.first].smutations,
                                         h.selected_keys, h.selected_row));
                    assert(validate_rows(pop.gametes[dip.second].smutations,
                                         h.selected_keys, h.selected_row2));
                    fill_matrix_with_rows(m, h, mtype);
                }
            // fill out other data fields
            for (auto &&i : neutral_keys)
                {
                    assert(pop.mutations[i.first].neutral);
                    m.neutral_positions.push_back(pop.mutations[i.first].pos);
                    m.neutral_popfreq.push_back(
                        static_cast<double>(pop.mcounts[i.first])
                        / static_cast<double>(2 * pop.diploids[deme].size()));
                }
            for (auto &&i : selected_keys)
                {
                    assert(!pop.mutations[i.first].neutral);
                    m.selected_positions.push_back(pop.mutations[i.first].pos);
                    m.selected_popfreq.push_back(
                        static_cast<double>(pop.mcounts[i.first])
                        / static_cast<double>(2 * pop.diploids[deme].size()));
                }
        }
        template <typename poptype>
        void
        fill_matrix(
            const poptype &pop, data_matrix &m,
            const std::vector<std::size_t> &individuals,
            const std::vector<std::pair<std::size_t, uint_t>> &neutral_keys,
            const std::vector<std::pair<std::size_t, uint_t>> &selected_keys,
            const std::size_t, sugar::MULTILOCPOP_TAG, matrix_type mtype)
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
                                              pop.gametes[locus.second], h);
                        }
                    fill_matrix_with_rows(m, h, mtype);
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
    }
}
