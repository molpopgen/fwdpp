#ifndef FWDPP_TS_GENERATE_DATA_MATRIX_HPP
#define FWDPP_TS_DATA_GENERATE_MATRIX_HPP

#include <vector>
#include <type_traits>
#include <cstdint>
#include <stdexcept>
#include "table_collection.hpp"
#include "marginal_tree_iterator.hpp"
#include <fwdpp/data_matrix.hpp>

namespace fwdpp
{
    namespace ts
    {
        namespace detail
        {
            inline void
            process_samples(const marginal_tree& marginal,
                            const TS_NODE_INT node, TS_NODE_INT index,
                            std::vector<std::int8_t>& genotypes)
            {
                auto right = marginal.right_sample[node];
                // Set all genotypes to ancestral state
                std::fill(genotypes.begin(), genotypes.end(), 0);
                while (true)
                    {
                        if (genotypes[index] == 1)
                            {
                                throw std::runtime_error("inconsist"
                                                         "ent "
                                                         "samples "
                                                         "list");
                            }
                        genotypes[index] = 1;
                        if (index == right)
                            {
                                break;
                            }
                        index = marginal.next_sample[index];
                    }
            }

            template <typename mcont_t>
            inline void
            update_data_matrix(const mcont_t& mutations, const std::size_t key,
                               const std::vector<std::int8_t>& genotypes,
                               data_matrix& rv)
            {
                auto n = mutations[key].neutral;
                auto& sm = (n) ? rv.neutral : rv.selected;
                auto& k = (n) ? rv.neutral_keys : rv.selected_keys;
                sm.positions.push_back(mutations[key].pos);
                sm.data.insert(sm.data.end(), genotypes.begin(),
                               genotypes.end());
                k.push_back(key);
            }
        } // namespace detail

        template <typename mcont_t>
        data_matrix
        generate_data_matrix(const table_collection& tables,
                             const std::vector<TS_NODE_INT>& samples,
                             const mcont_t& mutations,
                             const bool record_neutral,
                             const bool record_selected)
        {
            auto mut = tables.mutation_table.cbegin();
            const auto mut_end = tables.mutation_table.cend();
            marginal_tree_iterator tv(tables, samples);
            std::vector<std::int8_t> genotypes(samples.size(), 0);
            data_matrix rv(samples.size());
            while (tv(std::true_type(), std::true_type()))
                {
                    // Advance the mutation table records until we are
                    // in the current tree
                    while (mut < mut_end
                           && mutations[mut->key].pos < tv.marginal.left)
                        {
                            ++mut;
                        }
                    // Process mutations on this tree
                    for (; mut < mut_end
                           && mutations[mut->key].pos < tv.marginal.right;
                         ++mut)
                        {
                            bool is_neutral = mutations[mut->key].neutral;
                            if ((is_neutral && record_neutral)
                                || (!is_neutral && record_selected))
                                {
                                    auto index
                                        = tv.marginal.left_sample[mut->node];
                                    // Check if mutation leads to a sample
                                    if (index != TS_NULL_NODE)
                                        {
                                            detail::process_samples(
                                                tv.marginal, mut->node, index,
                                                genotypes);
                                            // Update our return value
                                            detail::update_data_matrix(
                                                mutations, mut->key, genotypes,
                                                rv);
                                        }
                                }
                        }
                    // Terminate early if we're lucky
                    if (!(mut < mut_end))
                        {
                            break;
                        }
                }
            return rv;
        }
    } // namespace ts
} // namespace fwdpp

#endif
