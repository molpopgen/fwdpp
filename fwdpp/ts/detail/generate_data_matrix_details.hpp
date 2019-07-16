#ifndef FWDPP_TS_GENERATE_DATA_MATRIX_DETAILS_HPP
#define FWDPP_TS_GENERATE_DATA_MATRIX_DETAILS_HPP

#include <cstdint>
#include <vector>
#include <algorithm>
#include <fwdpp/data_matrix.hpp>
#include "../marginal_tree.hpp"
#include "../marginal_tree_functions/samples.hpp"
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
                bool record_neutral, record_selected, skip_fixed;
                std::vector<std::int8_t> genotypes;
                data_matrix dm;

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
                data_matrix_filler(std::size_t samplesize, bool neutral,
                                   bool selected, bool nofixed)
                    : record_neutral(neutral), record_selected(selected),
                      skip_fixed(nofixed), genotypes(samplesize, -1),
                      dm(samplesize)
                {
                }

                inline void
                operator()(const marginal_tree& tree, const site& current_site,
                           mutation_key_vector::const_iterator mut,
                           mutation_key_vector::const_iterator mut_end)
                {
                    std::fill(genotypes.begin(), genotypes.end(),
                              current_site.ancestral_state);
                    int neutral = -1, selected = -1;
                    int nsamples = 0;
                    for (; mut < mut_end; ++mut)
                        {
                            neutral += (mut->neutral == true);
                            selected += (mut->neutral == false);
                            std::size_t lc = tree.leaf_counts[mut->node];
                            if ((mut->neutral && record_neutral)
                                || (selected && record_selected))
                                {
                                    if (lc > 0
                                        && (!skip_fixed || (lc < dm.ncol)))
                                        {
                                            process_samples(
                                                tree,
                                                convert_sample_index_to_nodes(
                                                    true),
                                                mut->node,
                                                [mut, &nsamples, this](
                                                    fwdpp::ts::TS_NODE_INT u) {
                                                    ++nsamples;
                                                    genotypes[u]
                                                        = mut->derived_state;
                                                });
                                        }
                                }
                        }
                    if (neutral != -1 && selected != -1)
                        {
                            throw tables_error("inconsistent neutral flags in "
                                               "mutation table");
                        }
                    if (nsamples)
                        {
                            update_data_matrix((mut_end - 1)->key,
                                               current_site.position,
                                               (neutral != -1) ? 1 : 0);
                        }
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
