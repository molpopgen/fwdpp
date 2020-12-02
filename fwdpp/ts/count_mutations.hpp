#ifndef FWDPP_TS_COUNT_MUTATIONS_HPP
#define FWDPP_TS_COUNT_MUTATIONS_HPP

#include <cassert>
#include <type_traits>
#include <vector>
#include "definitions.hpp"
#include "tree_visitor.hpp"

namespace fwdpp
{
    namespace ts
    {
        template <typename TableCollectionType, typename SAMPLES,
                  typename MutationContainerType>
        void
        count_mutations(const TableCollectionType& tables,
                        const MutationContainerType& mutations, SAMPLES&& samples,
                        std::vector<std::uint32_t>& mcounts)
        {
            // Use Kelleher et al. (2016)'s Algorithm L
            // to march through each marginal tree and its leaf
            // counts. At the same time, we march through our mutation
            // table, which is sorted by position.
            std::fill(mcounts.begin(), mcounts.end(), 0);
            mcounts.resize(mutations.size(), 0);

            auto mtable_itr = tables.mutations.begin();
            auto mtable_end = tables.mutations.end();
            tree_visitor<TableCollectionType> mti(tables, std::forward<SAMPLES>(samples),
                                                  update_samples_list(false));
            while (mti())
                {
                    auto& tree = mti.tree();
                    while (mtable_itr < mtable_end
                           && mutations[mtable_itr->key].pos < tree.left)
                        {
                            ++mtable_itr;
                        }
                    while (mtable_itr < mtable_end
                           && mutations[mtable_itr->key].pos < tree.right)
                        {
                            assert(mutations[mtable_itr->key].pos >= tree.left);
                            assert(mutations[mtable_itr->key].pos < tree.right);
                            mcounts[mtable_itr->key]
                                = tree.leaf_counts[mtable_itr->node];
                            ++mtable_itr;
                        }
                }
        }

        template <typename TableCollectionType, typename SAMPLES,
                  typename MutationContainerType>
        void
        count_mutations(const TableCollectionType& tables,
                        const MutationContainerType& mutations, SAMPLES&& samples,
                        SAMPLES&& preserved_nodes, std::vector<std::uint32_t>& mcounts,
                        std::vector<std::uint32_t>& acounts)
        {
            // Use Kelleher et al. (2016)'s Algorithm L
            // to march through each marginal tree and its leaf
            // counts. At the same time, we march through our mutation
            // table, which is sorted by position.
            std::fill(mcounts.begin(), mcounts.end(), 0);
            mcounts.resize(mutations.size(), 0);
            std::fill(acounts.begin(), acounts.end(), 0);
            acounts.resize(mutations.size(), 0);

            auto mtable_itr = tables.mutations.begin();
            auto mtable_end = tables.mutations.end();
            tree_visitor<TableCollectionType> mti(tables, std::forward<SAMPLES>(samples),
                                                  std::forward<SAMPLES>(preserved_nodes),
                                                  update_samples_list(false));
            while (mti())
                {
                    auto& tree = mti.tree();
                    while (mtable_itr < mtable_end
                           && mutations[mtable_itr->key].pos < tree.left)
                        {
                            ++mtable_itr;
                        }
                    while (mtable_itr < mtable_end
                           && mutations[mtable_itr->key].pos < tree.right)
                        {
                            assert(mutations[mtable_itr->key].pos >= tree.left);
                            assert(mutations[mtable_itr->key].pos < tree.right);
                            mcounts[mtable_itr->key]
                                = tree.leaf_counts[mtable_itr->node];
                            acounts[mtable_itr->key]
                                = tree.preserved_leaf_counts[mtable_itr->node];
                            ++mtable_itr;
                        }
                }
        }
    } // namespace ts
} // namespace fwdpp

#endif
