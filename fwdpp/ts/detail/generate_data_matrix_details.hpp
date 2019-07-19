#ifndef FWDPP_TS_GENERATE_DATA_MATRIX_DETAILS_HPP
#define FWDPP_TS_GENERATE_DATA_MATRIX_DETAILS_HPP

#include <cstdint>
#include <vector>
#include <algorithm>
#include <fwdpp/data_matrix.hpp>
#include "../marginal_tree.hpp"
#include "../site_visitor.hpp"
#include "../marginal_tree_functions/samples.hpp"
#include "../table_types.hpp"

namespace fwdpp
{
    namespace ts
    {
        namespace detail
        {
            template<typename SITE_CONST_ITER, typename MUT_CONST_ITR>
            inline void
            process_site_range(const marginal_tree & tree,
                               const SITE_CONST_ITER current_site,
                               const std::pair<MUT_CONST_ITR,MUT_CONST_ITR> & muts,
                               bool record_neutral, bool record_selected,
                               bool skip_fixed,
                               std::vector<std::int8_t>& genotypes,
                               data_matrix& dm)
            {
                int neutral = -1, selected = -1;
                std::fill(begin(genotypes), end(genotypes),
                          current_site->ancestral_state);
                int nsamples = 0;
                convert_sample_index_to_nodes convert(false);
                for (auto mut = muts.first; mut < muts.second; ++mut)
                    {
                        neutral += (mut->neutral == true);
                        selected += (mut->neutral == false);
                        std::size_t lc = tree.leaf_counts[mut->node];

                        if ((mut->neutral && record_neutral)
                            || (!mut->neutral && record_selected))
                            {
                                if (lc > 0
                                    && (!skip_fixed
                                        || (lc < genotypes.size())))
                                    {
                                        const auto f
                                            = [mut, &nsamples, &genotypes](
                                                  fwdpp::ts::TS_NODE_INT u) {
                                                  ++nsamples;
                                                  genotypes[u]
                                                      = mut->derived_state;
                                              };
                                        process_samples(tree, convert,
                                                        mut->node, f);
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
                        if (neutral != -1)
                            {
                                dm.neutral.positions.push_back(
                                    current_site->position);
                                dm.neutral_keys.push_back(
                                    (muts.second - 1)->key);
                                dm.neutral.data.insert(end(dm.neutral.data),
                                                       begin(genotypes),
                                                       end(genotypes));
                            }
                        else
                            {
                                dm.selected.positions.push_back(
                                    current_site->position);
                                dm.selected_keys.push_back(
                                    (muts.second - 1)->key);
                                dm.selected.data.insert(end(dm.selected.data),
                                                        begin(genotypes),
                                                        end(genotypes));
                            }
                    }
            }
        } // namespace detail
    }     // namespace ts
} // namespace fwdpp

#endif
