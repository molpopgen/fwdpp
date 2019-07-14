#ifndef FWDPP_TS_GENERATE_DATA_MATRIX_DETAILS_HPP
#define FWDPP_TS_GENERATE_DATA_MATRIX_DETAILS_HPP

#include <cstdint>
#include <vector>
#include <algorithm>
#include <fwdpp/data_matrix.hpp>
#include "../marginal_tree.hpp"
#include "../table_types.hpp"

namespace fwdpp
{
    namespace ts
    {
        namespace detail
        {
            class data_matrix_filler
            {
              private:
                site_vector::const_iterator site_table_begin;
                bool record_neutral, record_selected, skip_fixed;
                std::vector<std::int8_t> genotypes;
                bool filled;
                data_matrix dm;

                inline void
                process_samples(const marginal_tree& marginal,
                                const TS_NODE_INT node, TS_NODE_INT index)
                {
                    auto right = marginal.right_sample[node];
                    // Set all genotypes to ancestral state
                    //std::fill(genotypes.begin(), genotypes.end(), 0);
                    int x = 0;
                    while (true)
                        {
                            ++x;
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

                inline void
                update_data_matrix(const std::size_t key,
                                   const double position, const bool neutral)
                {
                    auto& sm = (neutral) ? dm.neutral : dm.selected;
                    auto& k = (neutral) ? dm.neutral_keys : dm.selected_keys;
                    sm.positions.push_back(position);
                    sm.data.insert(sm.data.end(), genotypes.begin(),
                                   genotypes.end());
                    k.push_back(key);
                }

              public:
                data_matrix_filler(site_vector::const_iterator b, bool rn,
                                   bool rs, bool sf, std::size_t sample_size)
                    : site_table_begin(b), record_neutral(rn),
                      record_selected(rs), skip_fixed(sf),
                      genotypes(sample_size, -1), filled(false),
                      dm(sample_size)
                {
                }

                inline mutation_key_vector::const_iterator
                operator()(const marginal_tree& tree, const site& current_site,
                           mutation_key_vector::const_iterator mut,
                           mutation_key_vector::const_iterator mut_end)
                {
                    bool initialized_site = false;
                    // mcopy prevents and off-by-one error in the return value
                    auto mcopy = mut;
                    for (; mut < mut_end
                           && (site_table_begin + mut->site)->position
                                  == current_site.position;
                         ++mut)
                        {
                            std::size_t total_leaf_count
                                = tree.preserved_leaf_counts[mut->node]
                                  + tree.leaf_counts[mut->node];
                            if (total_leaf_count
                                && (!skip_fixed
                                    || (skip_fixed
                                        && total_leaf_count
                                               < tree.sample_size())))
                                {
                                    if ((mut->neutral && record_neutral)
                                        || (!mut->neutral && record_selected))
                                        {
                                            if (!initialized_site)
                                                {
                                                    std::fill(
                                                        begin(genotypes),
                                                        end(genotypes),
                                                        current_site
                                                            .ancestral_state);
                                                }
                                            auto index
                                                = tree.left_sample[mut->node];
                                            // Check if mutation leads to a sample
                                            if (index != TS_NULL_NODE)
                                                {
                                                    process_samples(tree,
                                                                    mut->node,
                                                                    index);
                                                    filled = true;
                                                }
                                            mcopy = mut;
                                        }
                                }
                            if (filled)
                                {
                                    update_data_matrix(mcopy->key,
                                                       current_site.position,
                                                       mcopy->neutral);
                                }
                        }
                    return mut;
                }

                data_matrix
                data()
                {
                    return data_matrix(std::move(dm));
                }
            };
        } // namespace detail
    }     // namespace ts
} // namespace fwdpp

#endif
