#ifndef FWDPP_TS_GENERATE_DATA_MATRIX_DETAILS_HPP
#define FWDPP_TS_GENERATE_DATA_MATRIX_DETAILS_HPP

#include <cstdint>
#include <vector>
#include <algorithm>
#include <fwdpp/data_matrix.hpp>
#include "../marginal_tree.hpp"
#include "../site.hpp"
#include "../mutation_record.hpp"

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
            update_data_matrix(const std::size_t key, const double position,
                               const bool neutral,
                               const std::vector<std::int8_t>& genotypes,
                               data_matrix& rv)
            {
                auto& sm = (neutral) ? rv.neutral : rv.selected;
                auto& k = (neutral) ? rv.neutral_keys : rv.selected_keys;
                sm.positions.push_back(position);
                sm.data.insert(sm.data.end(), genotypes.begin(),
                               genotypes.end());
                k.push_back(key);
            }

            template <typename site_table_iterator, typename mtable_iterator>
            std::pair<bool, mtable_iterator>
            process_site(const marginal_tree& tree, const site& current_site,
                         site_table_iterator site_table_begin,
                         mtable_iterator mut, mtable_iterator mut_end,
                         bool record_neutral, bool record_selected,
                         bool skip_fixed, std::vector<std::int8_t>& genotypes)
            {
                bool initialized_site = false;
                bool filled = false;
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
                                    && total_leaf_count < tree.sample_size())))
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
                                                process_samples(
                                                    tree, mut->node, index,
                                                    genotypes);
                                                filled = true;
                                            }
                                    }
                            }
                    }
                return std::make_pair(filled, mut);
            }
        } // namespace detail
    }     // namespace ts
} // namespace fwdpp

#endif
