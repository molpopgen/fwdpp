#ifndef FWDPP_TS_TABLE_COLLECTION_HPP
#define FWDPP_TS_TABLE_COLLECTION_HPP

#include <vector>
#include <utility>
#include <cstdint>
#include <cmath>
#include <tuple>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <stdexcept>
#include <fwdpp/forward_types.hpp>
#include "node.hpp"
#include "edge.hpp"
#include "msprime_algo.hpp"
#include "mutation_record.hpp"
#include "indexed_edge.hpp" //TODO: create fewer header dependencies

namespace fwdpp
{
    namespace ts
    {
        using edge_vector = std::vector<edge>;
        using node_vector = std::vector<node>;
        using mutation_key_vector = std::vector<mutation_record>;
        struct table_collection
        {
          private:
            edge_vector temp_edges; //used for sorting
            void
            split_breakpoints_add_edges(
                const std::vector<double>& breakpoints,
                const std::tuple<TS_NODE_INT, TS_NODE_INT>& parents,
                const TS_NODE_INT next_index)
            {
                if (breakpoints.front() != 0.0)
                    {
                        this->push_back_edge(0., breakpoints.front(),
                                             std::get<0>(parents), next_index);
                    }
                for (unsigned j = 1; j < breakpoints.size(); ++j)
                    {
                        double a = breakpoints[j - 1];
                        double b = (j < breakpoints.size() - 1)
                                       ? breakpoints[j]
                                       : L;
                        if (b <= a)
                            {
                                throw std::runtime_error(
                                    "right must be > left");
                            }
                        if (j % 2 == 0.)
                            {
                                this->push_back_edge(
                                    a, b, std::get<0>(parents), next_index);
                            }
                        else
                            {
                                this->push_back_edge(
                                    a, b, std::get<1>(parents), next_index);
                            }
                    }
            }

            void
            split_breakpoints(
                const std::vector<double>& breakpoints,
                const std::tuple<TS_NODE_INT, TS_NODE_INT>& parents,
                const TS_NODE_INT next_index)
            {
                if (breakpoints.empty())
                    {
                        this->push_back_edge(0., L, std::get<0>(parents),
                                             next_index);
                        return;
                    }
                auto itr = std::adjacent_find(std::begin(breakpoints),
                                              std::end(breakpoints));
                if (itr == std::end(breakpoints))
                    {
                        split_breakpoints_add_edges(breakpoints, parents,
                                                    next_index);
                    }
                else
                    {
                        // Here, we need to reduce the input
                        // breakpoints to only those seen
                        // an odd number of times.
                        // Even numbers of the same breakpoint
                        // are "double x-overs" and thus
                        // cannot affect the genealogy.
                        std::vector<double> odd_breakpoints;
                        auto start = breakpoints.begin();
                        while (itr < breakpoints.end())
                            {
                                auto not_equal
                                    = std::find_if(itr, breakpoints.end(),
                                                   [itr](const double d) {
                                                       return d != *itr;
                                                   });
                                int even = (std::distance(itr, not_equal) % 2
                                            == 0.0);
                                odd_breakpoints.insert(odd_breakpoints.end(),
                                                       start, itr + 1 - even);
                                start = not_equal;
                                itr = std::adjacent_find(
                                    start, std::end(breakpoints));
                            }
                        odd_breakpoints.insert(odd_breakpoints.end(), start,
                                               breakpoints.end());
                        split_breakpoints_add_edges(odd_breakpoints, parents,
                                                    next_index);
                    }
            }

          public:
            node_vector node_table;
            edge_vector edge_table;
            mutation_key_vector mutation_table;
            indexed_edge_container input_left, output_right;
            // This reflects the length of
            // tables.edge_table after last simplification.
            // It can be used to make sure we only sort
            // newly-added nodes.
            std::ptrdiff_t edge_offset;
            const double L;
            std::vector<TS_NODE_INT> preserved_nodes;
            table_collection(const double maxpos)
                : temp_edges{}, node_table{}, edge_table{}, mutation_table{},
                  input_left{}, output_right{}, edge_offset{ 0 }, L{ maxpos },
                  preserved_nodes{}
            {
                if (maxpos < 0 || !std::isfinite(maxpos))
                    {
                        throw std::invalid_argument(
                            "maxpos must be > 0 and finite");
                    }
            }

            table_collection(const TS_NODE_INT num_initial_nodes,
                             const double initial_time, TS_NODE_INT pop,
                             const double maxpos)
                : temp_edges{}, node_table{}, edge_table{}, mutation_table{},
                  input_left{}, output_right{}, edge_offset{ 0 }, L{ maxpos },
                  preserved_nodes{}
            {
                if (maxpos < 0 || !std::isfinite(maxpos))
                    {
                        throw std::invalid_argument(
                            "maxpos must be > 0 and finite");
                    }
                for (TS_NODE_INT i = 0; i < num_initial_nodes; ++i)
                    {
                        node_table.push_back(node{ pop, initial_time });
                    }
            }

            void
            sort_edges() noexcept
            /// Sort the edge table.  On PARENT birth times.
            /// The sorting differs from msprime here. The difference
            /// is that we  assume that birth times are recorded forward in
            /// time rather than backwards.
            {
                std::sort(edge_table.begin() + edge_offset, edge_table.end(),
                          [this](const edge& a, const edge& b) {
                              auto ga = this->node_table[a.parent].generation;
                              auto gb = this->node_table[b.parent].generation;
                              if (ga == gb)
                                  {
                                      if (a.parent == b.parent)
                                          {
                                              if (a.child == b.child)
                                                  {
                                                      return a.left < b.left;
                                                  }
                                              return a.child < b.child;
                                          }
                                      return a.parent < b.parent;
                                  }
                              return ga > gb;
                          });
                if (edge_offset > 0)
                    {
                        temp_edges.reserve(edge_table.size());
#ifndef NDEBUG
                        auto size = edge_table.size();
#endif
                        temp_edges.clear();
                        temp_edges.insert(
                            temp_edges.end(),
                            std::make_move_iterator(edge_table.begin()
                                                    + edge_offset),
                            std::make_move_iterator(edge_table.end()));
                        temp_edges.insert(
                            temp_edges.end(),
                            std::make_move_iterator(edge_table.begin()),
                            std::make_move_iterator(edge_table.begin()
                                                    + edge_offset));
                        assert(temp_edges.size() == size);
                        temp_edges.swap(edge_table);
                    }
                temp_edges.clear();
                // TODO: allow for exceptions
                // rather than assertions.
                assert(edges_are_sorted());
            }

            template <typename mutation_container>
            void
            sort_mutations(const mutation_container& mutations)
            {
                //mutations are sorted by increasing position
                std::sort(mutation_table.begin(), mutation_table.end(),
                          [&mutations](const mutation_record& a,
                                       const mutation_record& b) {
                              return mutations[a.key].pos
                                     < mutations[b.key].pos;
                          });
            }

            template <typename mutation_container>
            void
            sort_tables(const mutation_container& mutations)
            /// Sorts the tables
            /// Note that mutations can be mocked via any struct
            /// containing double pos
            {
                sort_edges();
                sort_mutations(mutations);
            }

            bool
            edges_are_sorted() const noexcept
            /// Test the MINIMAL sorting requirement.
            /// This minimal condition is important,
            /// as ancient sample tracking will not
            /// guarantee a sort order on the parental
            /// node IDs.
            {
                return std::is_sorted(
                    edge_table.begin(), edge_table.end(),
                    [this](const edge& a, const edge& b) {
                        auto ga = this->node_table[a.parent].generation;
                        auto gb = this->node_table[b.parent].generation;
                        return ga > gb
                               && (std::tie(a.child, a.left)
                                   < std::tie(b.child, b.left));
                    });
            }

            void
            clear() noexcept
            /// Clears internal vectors.
            /// Mostly used during simplification
            /// where a table_collection is
            /// used as a temp object.
            {
                node_table.clear();
                edge_table.clear();
                mutation_table.clear();
                preserved_nodes.clear();
            }

            void
            record_preserved_nodes(const std::vector<TS_NODE_INT>& node_ids)
            {
                for (auto i : node_ids)
                    {
                        if (i >= node_table.size())
                            {
                                throw std::invalid_argument(
                                    "node id larger than node table size");
                            }
                        if (std::find(preserved_nodes.begin(),
                                      preserved_nodes.end(), i)
                            != preserved_nodes.end())
                            {
                                throw std::invalid_argument(
                                    "node already recorded as a "
                                    "preserved_node");
                            }
                        preserved_nodes.push_back(i);
                    }
            }

            void
            push_back_node(double generation, std::int32_t pop)
            {
                node_table.push_back(node{ pop, generation });
            }

            template <typename... args>
            void
            emplace_back_node(args&&... Args)
            {
                node_table.emplace_back(node{ std::forward<args>(Args)... });
            }

            void
            push_back_edge(double l, double r, TS_NODE_INT parent,
                           TS_NODE_INT child)
            {
                edge_table.push_back(edge{ l, r, parent, child });
            }

            template <typename... args>
            void
            emplace_back_edge(args&&... Args)
            {
                edge_table.emplace_back(edge{ std::forward<args>(Args)... });
            }

            void
            build_indexes()
            /// Generates the index vectors referred to
            /// as I and O in Kelleher et al. (2016)
            {
                input_left.reserve(edge_table.size());
                output_right.reserve(edge_table.size());
                input_left.clear();
                output_right.clear();
                for (auto& e : edge_table)
                    {
                        assert(e.left < e.right);
                        input_left.emplace_back(
                            e.left, -node_table[e.parent].generation, e.parent,
                            e.child);
                        output_right.emplace_back(
                            e.right, node_table[e.parent].generation, e.parent,
                            e.child);
                    }
                std::sort(input_left.begin(), input_left.end());
                std::sort(output_right.begin(), output_right.end());
            }

            void
            add_offspring_data(
                const TS_NODE_INT next_index,
                const std::vector<double>& breakpoints,
                const std::tuple<TS_NODE_INT, TS_NODE_INT>& parents,
                const double generation)
            {
                emplace_back_node(0, generation);
                split_breakpoints(breakpoints, parents, next_index);
            }

            void
            add_offspring_data(
                const TS_NODE_INT next_index,
                const std::vector<double>& breakpoints,
                const std::vector<std::uint32_t>& new_mutations,
                const std::tuple<TS_NODE_INT, TS_NODE_INT>& parents,
                const double generation)
            {
                add_offspring_data(next_index, breakpoints, parents,
                                   generation);
                for (auto& m : new_mutations)
                    {
                        mutation_table.emplace_back(
                            mutation_record{ next_index, m });
                        assert(mutation_table.back().node == next_index);
                        assert(mutation_table.back().key == m);
                    }
            }

            std::size_t
            num_nodes() const
            {
                return node_table.size();
            }

            void
            update_offset()
            {
                edge_offset = edge_table.size();
            }

            template <typename mcont_t>
            void
            count_mutations(const mcont_t& mutations,
                            const std::vector<TS_NODE_INT>& samples,
                            std::vector<fwdpp::uint_t>& mcounts)
            {
                // Use Kelleher et al. (2016)'s Algorithm L
                // to march through each marginal tree and its leaf
                // counts. At the same time, we march through our mutation
                // table, which is sorted by position.
                std::fill(mcounts.begin(), mcounts.end(), 0);
                mcounts.resize(mutations.size(), 0);

                auto mtable_itr = mutation_table.begin();
                auto mtable_end = mutation_table.end();
                auto mutation_counter = [&mutations, &mtable_itr, mtable_end,
                                         &mcounts](
                                            const marginal_tree& marginal) {
                    while (mtable_itr < mtable_end
                           && mutations[mtable_itr->key].pos < marginal.left)
                        {
                            ++mtable_itr;
                        }
                    while (mtable_itr < mtable_end
                           && mutations[mtable_itr->key].pos < marginal.right)
                        {
                            assert(mutations[mtable_itr->key].pos
                                   >= marginal.left);
                            assert(mutations[mtable_itr->key].pos
                                   < marginal.right);
                            mcounts[mtable_itr->key]
                                = marginal.leaf_counts[mtable_itr->node];
                            ++mtable_itr;
                        }
                };

                algorithmL(input_left, output_right, samples,
                           node_table.size(), L, mutation_counter);
            }
            template <typename mcont_t>
            void
            count_mutations(const mcont_t& mutations,
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

                auto mtable_itr = mutation_table.begin();
                auto mtable_end = mutation_table.end();
                auto mutation_counter = [&mutations, &mtable_itr, mtable_end,
                                         &mcounts,&acounts](
                                            const marginal_tree& marginal) {
                    while (mtable_itr < mtable_end
                           && mutations[mtable_itr->key].pos < marginal.left)
                        {
                            ++mtable_itr;
                        }
                    while (mtable_itr < mtable_end
                           && mutations[mtable_itr->key].pos < marginal.right)
                        {
                            assert(mutations[mtable_itr->key].pos
                                   >= marginal.left);
                            assert(mutations[mtable_itr->key].pos
                                   < marginal.right);
                            mcounts[mtable_itr->key]
                                = marginal.leaf_counts[mtable_itr->node];
                            acounts[mtable_itr->key]
                                = marginal.preserved_leaf_counts[mtable_itr->node];
                            ++mtable_itr;
                        }
                };

                algorithmL(input_left, output_right, samples, preserved_nodes,
                           node_table.size(), L, mutation_counter);
            }
        };
    } // namespace ts
} // namespace fwdpp

#endif
