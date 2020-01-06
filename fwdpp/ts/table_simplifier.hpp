#ifndef FWDPP_TS_TABLE_SIMPLIFIER_HPP
#define FWDPP_TS_TABLE_SIMPLIFIER_HPP

#include <iostream>
#include <cmath>
#include <map>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cstddef>
#include <stdexcept>
#include <unordered_map>
#include "node.hpp"
#include "edge.hpp"
#include "table_collection.hpp"

namespace fwdpp
{
    namespace ts
    {
        class table_simplifier
        /*! \brief Implements the simplification algorithm of \cite Kelleher2018-fu
         *
         *  Keeping a persistent simplifier around during a simplification avoids many 
         *  repeated large memory allocations.  On a single-core machine, the more often
         *  you simplify, the more this matters.  When running many simulations on a many-core 
         *  machine, keeping a persistent object reduces competition between threads for big
         *  memory chunks, which should also lead to a performance boost.
         *
         *  Many of the implementation details are private functions, which are subject to change
         *  without notice.
         *
         *  \version 0.7.0 Added to fwdpp
         */
        {
          private:
            struct segment
            {
                double left, right;
                TS_NODE_INT node;
                segment(double l, double r, TS_NODE_INT n)
                    : left{ l }, right{ r }, node{ n }
                {
                    if (right <= left)
                        {
                            throw std::invalid_argument(
                                "right must be > left");
                        }
                }
            };

            struct mutation_node_map_entry
            {
                TS_NODE_INT node;
                std::size_t site, location;
                mutation_node_map_entry(TS_NODE_INT n, std::size_t s,
                                        std::size_t l)
                    : node(n), site(s), location(l)
                {
                }
            };

            class segment_overlapper
            /// This class is an iterable object
            /// over [left, right) -> segment
            /// mappings, where the segments
            /// are the genomic intervals in
            /// child nodes overlapping with
            /// the current parent in the
            /// current genomic interval.
            {
              private:
                std::vector<segment>::const_iterator sbeg, send;

                inline double
                set_partition()
                {
                    double tright = std::numeric_limits<double>::max();
                    auto b = overlapping.begin();
                    for (auto i = overlapping.begin(); i < overlapping_end;
                         ++i)
                        {
                            if (i->right > left)
                                {
                                    *b = *i;
                                    tright = std::min(tright, b->right);
                                    ++b;
                                }
                        }
                    overlapping_end = b;
                    return tright;
                }

              public:
                std::vector<segment> overlapping;
                std::vector<segment>::iterator overlapping_end;
                double left, right;
                segment_overlapper()
                    : sbeg(), send(), overlapping{},
                      overlapping_end(overlapping.end()), left(0),
                      right(std::numeric_limits<double>::max())
                {
                }

                void
                init(std::vector<segment>& segs)
                {
                    sbeg = segs.begin();
                    // The - 1 for send assumes a "cap"/sentinel value.
                    send = segs.end() - 1;
                    overlapping.clear();
                    overlapping_end = overlapping.end();
                    left = 0.0;
                    right = std::numeric_limits<double>::max();
                }

                bool
                operator()()
                {
                    bool rv = 0;
                    if (sbeg < send)
                        {
                            left = right;
                            auto tright = set_partition();
                            if (num_overlaps() == 0)
                                {
                                    left = sbeg->left;
                                }
                            while (sbeg < send && sbeg->left == left)
                                {
                                    tright = std::min(tright, sbeg->right);
                                    overlapping_end
                                        = overlapping.insert(overlapping_end,
                                                             *sbeg)
                                          + 1;
                                    ++sbeg;
                                }
                            right = std::min(sbeg->left, tright);
                            rv = true;
                        }
                    else
                        {
                            left = right;
                            right = std::numeric_limits<double>::max();
                            auto tright = set_partition();
                            if (num_overlaps() > 0)
                                {
                                    right = tright;
                                    rv = true;
                                }
                        }
                    return rv;
                }

                std::int64_t
                num_overlaps()
                {
                    return std::distance(overlapping.begin(), overlapping_end);
                }
            };
            // These are temp tables/buffer
            // for simplification.  We keep
            // their allocated memory persistent.
            edge_vector new_edge_table;
            node_vector new_node_table;
            site_vector new_site_table;
            // segment_queue mimics a min queue of segments w.r.to
            // segment::left.
            std::vector<segment> segment_queue;
            std::vector<std::vector<segment>> Ancestry;
            /// Temp container used for compacting edges
            edge_vector E;
            // region length
            const double L;
            segment_overlapper o;
            std::vector<mutation_node_map_entry> mutation_map;

            void
            cleanup() noexcept
            // Clears out data from
            // temp containers after simplify.
            // Retains container capacity.
            {
                new_edge_table.clear();
                new_node_table.clear();
                new_site_table.clear();
                E.clear();
                // It is tempting to
                // just clear out each
                // inner element. Bad idea.
                // You use >= 10X more RAM
                // in "big" simulations.
                Ancestry.clear();
            }

            edge_vector::const_iterator
            step_S3(edge_vector::const_iterator edge_ptr,
                    const edge_vector::const_iterator edge_end, TS_NODE_INT u)
            {
                segment_queue.clear();
                for (; edge_ptr < edge_end && edge_ptr->parent == u;
                     ++edge_ptr)
                    {
                        // For each edge corresponding to this parent,
                        // we look at all segments from the child.
                        // If the two segments overlap, we add the
                        // minimal
                        // overlap to our queue.
                        // This is Step S3.
                        for (auto& seg : Ancestry[edge_ptr->child])
                            {
                                if (seg.right > edge_ptr->left
                                    && edge_ptr->right > seg.left)
                                    {
                                        segment_queue.emplace_back(
                                            std::max(seg.left, edge_ptr->left),
                                            std::min(seg.right,
                                                     edge_ptr->right),
                                            seg.node);
                                    }
                            }
                    }
                // Sort for processing via the overlapper
                std::sort(segment_queue.begin(), segment_queue.end(),
                          [](const segment& a, const segment& b) {
                              return a.left < b.left;
                          });
                // Add sentinel
                segment_queue.emplace_back(
                    segment{ L, L + 1.0, TS_NULL_NODE });
                return edge_ptr;
            }

            void
            buffer_edge(std::vector<edge>& buffered_edges, const double left,
                        const double right, const TS_NODE_INT parent,
                        const TS_NODE_INT child)
            {
                auto itr = std::find_if(
                    buffered_edges.rbegin(), buffered_edges.rend(),
                    [child](const edge& e) { return e.child == child; });
                if (itr == buffered_edges.rend())
                    {
                        buffered_edges.emplace_back(
                            edge{ left, right, parent, child });
                    }
                else
                    {
                        if (itr->right == left)
                            {
                                itr->right = right;
                            }
                        else
                            {
                                buffered_edges.emplace_back(
                                    edge{ left, right, parent, child });
                            }
                    }
            }

            void
            add_ancestry(TS_NODE_INT input_id, double left, double right,
                         TS_NODE_INT node)
            {
                if (Ancestry[input_id].empty())
                    {
                        Ancestry[input_id].emplace_back(left, right, node);
                    }
                else
                    {
                        auto& last = Ancestry[input_id].back();
                        if (last.right == left && last.node == node)
                            {
                                last.right = right;
                            }
                        else
                            {
                                Ancestry[input_id].emplace_back(left, right,
                                                                node);
                            }
                    }
            }

            void
            merge_ancestors(const node_vector& input_node_table,
                            const TS_NODE_INT parent_input_id,
                            std::vector<TS_NODE_INT>& idmap)
            {
                auto output_id = idmap[parent_input_id];
                bool is_sample = (output_id != TS_NULL_NODE);
                if (is_sample == true)
                    {
                        Ancestry[parent_input_id].clear();
                    }
                double previous_right = 0.0;
                o.init(segment_queue);
                TS_NODE_INT ancestry_node = TS_NULL_NODE;
                E.clear();
                while (o() == true)
                    {
                        if (o.num_overlaps() == 1)
                            {
                                ancestry_node = o.overlapping[0].node;
                                if (is_sample)
                                    {
                                        buffer_edge(E, o.left, o.right,
                                                    output_id, ancestry_node);
                                        ancestry_node = output_id;
                                    }
                            }
                        else
                            {
                                if (output_id == TS_NULL_NODE)
                                    {
                                        new_node_table.emplace_back(node{
                                            input_node_table[parent_input_id]
                                                .deme,
                                            input_node_table[parent_input_id]
                                                .time });
                                        output_id = new_node_table.size() - 1;
                                        // update sample map
                                        idmap[parent_input_id] = output_id;
                                    }
                                ancestry_node = output_id;
                                for (auto x = o.overlapping.begin();
                                     x < o.overlapping_end; ++x)
                                    {
                                        buffer_edge(E, o.left, o.right,
                                                    output_id, x->node);
                                    }
                            }
                        if (is_sample && o.left != previous_right)
                            {
                                add_ancestry(parent_input_id, previous_right,
                                             o.left, output_id);
                            }
                        add_ancestry(parent_input_id, o.left, o.right,
                                     ancestry_node);
                        previous_right = o.right;
                    }
                if (is_sample && previous_right != L)
                    {
                        add_ancestry(parent_input_id, previous_right, L,
                                     output_id);
                    }
                if (output_id != TS_NULL_NODE)
                    {
                        auto n = output_buffered_edges(E);
                        if (!n && !is_sample)
                            {
                                new_node_table.erase(new_node_table.begin()
                                                         + output_id,
                                                     new_node_table.end());
                                idmap[parent_input_id] = TS_NULL_NODE;
                            }
                    }
            }

            std::size_t
            output_buffered_edges(std::vector<edge>& buffered_edges)
            /// Take our buffered edges and add them to the output edge table
            {
                std::stable_sort(buffered_edges.begin(), buffered_edges.end(),
                                 [](const edge& a, const edge& b) {
                                     return a.child < b.child;
                                 });
                new_edge_table.insert(new_edge_table.end(),
                                      buffered_edges.begin(),
                                      buffered_edges.end());
                return buffered_edges.size();
            }

            void
            prep_mutation_simplification(
                const site_vector& sites,
                const mutation_key_vector& mutation_table)
            {
                mutation_map.clear();
                mutation_map.reserve(mutation_table.size());
                for (std::size_t i = 0; i < mutation_table.size(); ++i)
                    {
                        mutation_map.emplace_back(mutation_table[i].node,
                                                  mutation_table[i].site, i);
                    }

                std::sort(mutation_map.begin(), mutation_map.end(),
                          [&sites](const mutation_node_map_entry& a,
                                   const mutation_node_map_entry& b) {
                              return std::tie(a.node, sites[a.site].position)
                                     < std::tie(b.node,
                                                sites[b.site].position);
                          });
            }

            void
            record_site(const site_vector& sites, mutation_record& mr)
            {
                double pos = sites[mr.site].position;
                if (new_site_table.empty()
                    || new_site_table.back().position != pos)
                    {
                        new_site_table.push_back(sites[mr.site]);
                    }
                mr.site = new_site_table.size() - 1;
            }

            std::vector<std::size_t>
            simplify_mutations(mutation_key_vector& mt, site_vector& sites)
            // Remove all mutations that do not map to nodes
            // in the simplified tree.  The key here is
            // that Ancestry contains the history of
            // each node, which we use for the remapping.
            {
                // Set all output nodes to null for now.
                for (auto& mr : mt)
                    {
                        mr.node = TS_NULL_NODE;
                    }

                // Map the input node id of a mutation to
                // its output node id.  If no output ID exists,
                // then the mutation will be removed by the
                // call to erase below.
                auto map_itr = mutation_map.begin();
                const auto map_end = mutation_map.end();

                while (map_itr < map_end)
                    {
                        auto n = map_itr->node;
                        auto seg = Ancestry[n].cbegin();
                        const auto seg_e = Ancestry[n].cend();
                        for (; map_itr < map_end
                               && map_itr->node == n;) //++map_itr)
                            {
                                if (seg == seg_e)
                                    {
                                        ++map_itr;
                                        break;
                                    }
                                while (seg < seg_e && map_itr < map_end
                                       && map_itr->node == n)
                                    {
                                        auto pos
                                            = sites[map_itr->site].position;
                                        if (seg->left <= pos
                                            && pos < seg->right)
                                            {
                                                mt[map_itr->location].node
                                                    = seg->node;
                                                ++map_itr;
                                            }
                                        else if (pos >= seg->right)
                                            {
                                                ++seg;
                                            }
                                        else
                                            {
                                                ++map_itr;
                                            }
                                    }
                            }
                    }

                // Any mutations with null node values do not have
                // ancestry and may be removed.
                std::vector<std::size_t> preserved_variants;
                auto itr = std::remove_if(mt.begin(), mt.end(),
                                          [](const mutation_record& mr) {
                                              return mr.node == TS_NULL_NODE;
                                          });
                preserved_variants.reserve(std::distance(itr, mt.end()));
                for (auto i = mt.begin(); i != itr; ++i)
                    {
                        record_site(sites, *i);
                        preserved_variants.push_back(i->key);
                    }

                mt.erase(itr, mt.end());
                sites.swap(new_site_table);
                //TODO: replace assert with exception
                assert(std::is_sorted(mt.begin(), mt.end(),
                                      [&sites](const mutation_record& a,
                                               const mutation_record& b) {
                                          return sites[a.site].position
                                                 < sites[b.site].position;
                                      }));
                return preserved_variants;
            }

            void
            record_sample_nodes(const std::vector<TS_NODE_INT>& samples,
                                const table_collection& tables,
                                std::vector<TS_NODE_INT>& idmap)
            /// \version 0.7.1 Throw exception if a sample is recorded twice
            {
                for (const auto& s : samples)
                    {
                        // See GitHub issue 158
                        // for background
                        if (idmap[s] != TS_NULL_NODE)
                            {
                                throw std::invalid_argument(
                                    "invalid sample list");
                            }
                        new_node_table.emplace_back(
                            node{ tables.node_table[s].deme,
                                  tables.node_table[s].time });
                        add_ancestry(s, 0, L,
                                     static_cast<TS_NODE_INT>(
                                         new_node_table.size() - 1));
                        idmap[s] = static_cast<TS_NODE_INT>(
                            new_node_table.size() - 1);
                    }
            }

          public:
            explicit table_simplifier(const double maxpos)
                : new_edge_table{}, new_node_table{}, new_site_table{},
                  segment_queue{}, Ancestry{}, E{}, L{ maxpos }, o{},
                  mutation_map{}
            {
                if (maxpos < 0 || !std::isfinite(maxpos))
                    {
                        throw std::invalid_argument(
                            "maxpos must be > 0 and finite");
                    }
            }

            std::pair<std::vector<TS_NODE_INT>, std::vector<std::size_t>>
            simplify(table_collection& tables,
                     const std::vector<TS_NODE_INT>& samples)
            /// Simplify algorithm is approximately the same
            /// logic as used in msprime 0.6.0
            ///
            /// \param tables A table_collection
            /// \param samples A list of sample (node) ids.
            /// \version 0.7.1 Throw exception if a sample is recorded twice
            /// \version 0.7.3 Return value is now a pair containing the
            /// node ID map and a vector of keys to mutations preserved in
            /// mutation tables
            {
                Ancestry.resize(tables.node_table.size());

                // Set some things up for later mutation simplification
                prep_mutation_simplification(tables.site_table,
                                             tables.mutation_table);

                // Relates input node ids to output node ids
                std::vector<TS_NODE_INT> idmap(tables.node_table.size(),
                                               TS_NULL_NODE);

                // We take our samples and add them to both the output
                // node list and initialize their ancestry with
                // a segment on [0,L).
                record_sample_nodes(samples, tables, idmap);
                // Add samples for any preserved nodes in the tables:
                record_sample_nodes(tables.preserved_nodes, tables, idmap);

                // At this point, our edges are sorted by birth
                // order of parents, from present to past.
                // We can now work our way up the pedigree.
                // This outer loop differs from how we describe it in the
                // paper, but the strict sorting of edges means that this
                // equivalent.
                auto edge_ptr = tables.edge_table.cbegin();
                const auto edge_end = tables.edge_table.cend();
                while (edge_ptr < edge_end)
                    {
                        auto u = edge_ptr->parent;
                        edge_ptr = step_S3(edge_ptr, edge_end, u);
                        merge_ancestors(tables.node_table, u, idmap);
                    }

                // When there are preserved nodes, we need to re map
                // their input ids to output ids
                for (auto& p : tables.preserved_nodes)
                    {
                        if (idmap[p] == TS_NULL_NODE)
                            {
                                throw std::runtime_error(
                                    "preserved node output id maps to null");
                            }
                        p = idmap[p];
                    }

                assert(
                    static_cast<std::size_t>(std::count_if(
                        idmap.begin(), idmap.end(),
                        [](const TS_NODE_INT i) { return i != TS_NULL_NODE; }))
                    == new_node_table.size());
                // Update the tables.  To keep memory use as sane as possible,
                // we use resize-and-move here.  In theory, we can also do
                // vector swaps, but that has a side-effect of keeping
                // far too much RAM allocated compared to what we need.
                tables.edge_table.resize(new_edge_table.size());
                std::move(new_edge_table.begin(), new_edge_table.end(),
                          tables.edge_table.begin());
                tables.node_table.resize(new_node_table.size());
                std::move(new_node_table.begin(), new_node_table.end(),
                          tables.node_table.begin());
                // TODO: allow for exception instead of assert
                assert(tables.edges_are_minimally_sorted());
                tables.update_offset();
                auto preserved_variants = simplify_mutations(
                    tables.mutation_table, tables.site_table);

                cleanup();
                std::pair<std::vector<TS_NODE_INT>, std::vector<std::size_t>>
                    rv(std::move(idmap), std::move(preserved_variants));
                return rv;
            }
        };
    } // namespace ts
} // namespace fwdpp

#endif
