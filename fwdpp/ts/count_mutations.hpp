#ifndef FWDPP_TS_COUNT_MUTATIONS_HPP
#define FWDPP_TS_COUNT_MUTATIONS_HPP

#include <cassert>
#include <type_traits>
#include <vector>
#include "definitions.hpp"
#include "table_collection.hpp"
#include "marginal_tree_iterator.hpp"

namespace fwdpp
{
    namespace ts
    {
        template <typename mcont_t>
        void
        count_mutations(const table_collection& tables,
                        const mcont_t& mutations,
                        const std::vector<TS_NODE_INT>& samples,
                        std::vector<fwdpp::uint_t>& mcounts)
        {
            // Use Kelleher et al. (2016)'s Algorithm L
            // to march through each marginal tree and its leaf
            // counts. At the same time, we march through our mutation
            // table, which is sorted by position.
            std::fill(mcounts.begin(), mcounts.end(), 0);
            mcounts.resize(mutations.size(), 0);

            auto mtable_itr = tables.mutation_table.begin();
            auto mtable_end = tables.mutation_table.end();
            marginal_tree_iterator mti(tables, samples);
            while (mti(std::true_type(), std::false_type()))
                {
                    while (mtable_itr < mtable_end
                           && mutations[mtable_itr->key].pos
                                  < mti.marginal.left)
                        {
                            ++mtable_itr;
                        }
                    while (mtable_itr < mtable_end
                           && mutations[mtable_itr->key].pos
                                  < mti.marginal.right)
                        {
                            assert(mutations[mtable_itr->key].pos
                                   >= mti.marginal.left);
                            assert(mutations[mtable_itr->key].pos
                                   < mti.marginal.right);
                            mcounts[mtable_itr->key]
                                = mti.marginal.leaf_counts[mtable_itr->node];
                            ++mtable_itr;
                        }
                }
        }

        template <typename mcont_t>
        void
        count_mutations(const table_collection& tables,
                        const mcont_t& mutations,
                        const std::vector<TS_NODE_INT>& samples,
                        std::vector<fwdpp::uint_t>& mcounts,
                        std::vector<fwdpp::uint_t>& acounts)
        {
            // Use Kelleher et al. (2016)'s Algorithm L
            // to march through each marginal tree and its leaf
            // counts. At the same time, we march through our mutation
            // table, which is sorted by position.
            std::fill(mcounts.begin(), mcounts.end(), 0);
            mcounts.resize(mutations.size(), 0);
            std::fill(acounts.begin(), acounts.end(), 0);
            acounts.resize(mutations.size(), 0);

            auto mtable_itr = tables.mutation_table.begin();
            auto mtable_end = tables.mutation_table.end();
            marginal_tree_iterator mti(tables, samples,
                                       tables.preserved_nodes);
            while (mti(std::true_type(), std::false_type()))
                {
                    while (mtable_itr < mtable_end
                           && mutations[mtable_itr->key].pos
                                  < mti.marginal.left)
                        {
                            ++mtable_itr;
                        }
                    while (mtable_itr < mtable_end
                           && mutations[mtable_itr->key].pos
                                  < mti.marginal.right)
                        {
                            assert(mutations[mtable_itr->key].pos
                                   >= mti.marginal.left);
                            assert(mutations[mtable_itr->key].pos
                                   < mti.marginal.right);
                            mcounts[mtable_itr->key]
                                = mti.marginal.leaf_counts[mtable_itr->node];
                            acounts[mtable_itr->key]
                                = mti.marginal
                                      .preserved_leaf_counts[mtable_itr->node];
                            ++mtable_itr;
                        }
                }
        }
    } // namespace ts
} // namespace fwdpp

#endif