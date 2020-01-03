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
#include <fwdpp/ts/exceptions.hpp>
#include "table_types.hpp"
#include "indexed_edge.hpp" //TODO: create fewer header dependencies

namespace fwdpp
{
    namespace ts
    {
        struct table_collection
        /*!
		 * \brief A collection of tables for a single simulation.
         *
         * \version 0.7.0 Added to fwdpp
		 */
        {
          private:
            /// Length of the genomic region.
            double L;
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
                // TODO: replace with exception via a debug mode
                assert(std::count(begin(breakpoints), end(breakpoints),
                                  std::numeric_limits<double>::max())
                       == 1);
                assert(breakpoints.back()
                       == std::numeric_limits<double>::max());
                for (unsigned j = 1; j < breakpoints.size(); ++j)
                    {
                        double a = breakpoints[j - 1];
                        double b = (j < breakpoints.size() - 1)
                                       ? breakpoints[j]
                                       : L;
                        if (b <= a)
                            {
                                throw std::invalid_argument(
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

            void
            record_site_during_rebuild(const site& s, mutation_record& mr)
            {
                if (site_table.empty()
                    || site_table.back().position != s.position)
                    {
                        site_table.push_back(s);
                    }
                mr.site = site_table.size() - 1;
            }

          public:
            /// Node table for this simulation
            node_vector node_table;
            /// Edge table for this simulation
            edge_vector edge_table;
            /// Mutation table for this simulation;
            mutation_key_vector mutation_table;
            /// Site table
            site_vector site_table;
            /// The input edge vector. "I" in \cite Kelleher2016-cb, page 13
            indexed_edge_container input_left;
            /// The output edge vector. "O" in \cite Kelleher2016-cb, page 13
            indexed_edge_container output_right;
            /// This reflects the length of
            /// tables.edge_table after last simplification.
            /// It can be used to make sure we only sort
            /// newly-added nodes.
            std::ptrdiff_t edge_offset;
            /// A vector of dead/ancient sample nodes
            std::vector<TS_NODE_INT> preserved_nodes;

            explicit table_collection(const double maxpos)
                : L{ maxpos }, node_table{}, edge_table{}, mutation_table{},
                  site_table{}, input_left{}, output_right{}, edge_offset{ 0 },
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
                : L{ maxpos }, node_table{}, edge_table{}, mutation_table{},
                  site_table{}, input_left{}, output_right{}, edge_offset{ 0 },
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
            ///
            /// In between simplifications, edges are appended to the table
            /// starting at edge_offset.  edges prior to that point are,
            /// by definition, sorted b/c they are the output of simplification.
            /// Thus, we only sort the new edges and then perform a left rotate
            /// w.r.to edge_offset to put the nodes in the correct final order.
            /// The rotation is exactly edge_table.size() swaps.
            {
                std::sort(edge_table.begin() + edge_offset, edge_table.end(),
                          [this](const edge& a, const edge& b) {
                              auto ga = this->node_table[a.parent].time;
                              auto gb = this->node_table[b.parent].time;
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
                        std::rotate(edge_table.begin(),
                                    edge_table.begin() + edge_offset,
                                    edge_table.end());
                    }
#ifndef NDEBUG
                // TODO: allow for exceptions
                // rather than assertions.
                if (edge_offset == 0)
                    {
                        assert(edges_are_sorted());
                    }
                else if (preserved_nodes.empty())
                    {
                        assert(edges_are_sorted());
                    }
                else
                    {
                        assert(edges_are_minimally_sorted());
                    }
#endif
            }

            void
            sort_mutations()
            {
                //mutations are sorted by increasing position
                std::sort(mutation_table.begin(), mutation_table.end(),
                          [this](const mutation_record& a,
                                 const mutation_record& b) {
                              return site_table[a.site].position
                                     < site_table[b.site].position;
                          });
            }

            void
            sort_tables_for_simplification()
            /// Sorts the tables for simplification, which means only
            /// sorting edge and mutation tables, as the site table
            /// will be rebuilt during simplification.
            {
                sort_edges();
                sort_mutations();
            }

            void
            sort_tables()
            /// Sort all tables.  The site table is rebuilt.
            {
                sort_edges();
                sort_mutations_rebuild_site_table();
            }

            bool
            edges_are_sorted() const noexcept
            /// Test that edges are in proper sort order for simiplification
            /// \version 0.8.0 Change from testing minimal requirement to requirement
            /// for simplification.
            {
                return std::is_sorted(
                    edge_table.begin(), edge_table.end(),
                    [this](const edge& a, const edge& b) {
                        auto ga = this->node_table[a.parent].time;
                        auto gb = this->node_table[b.parent].time;
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
            }

            bool
            edges_are_minimally_sorted() const noexcept
            /// Test the MINIMAL sorting requirement.
            /// This minimal condition is important,
            /// as ancient sample tracking will not
            /// guarantee a sort order on the parental
            /// node IDs.
            /// \version 0.8.0 Added to library
            {
                return std::is_sorted(
                    edge_table.begin(), edge_table.end(),
                    [this](const edge& a, const edge& b) {
                        auto ga = this->node_table[a.parent].time;
                        auto gb = this->node_table[b.parent].time;
                        return ga > gb
                               && (std::tie(a.child, a.left)
                                   < std::tie(b.child, b.left));
                    });
            }

            void
            clear() noexcept
            /// Clears internal vectors.
            {
                node_table.clear();
                edge_table.clear();
                mutation_table.clear();
                preserved_nodes.clear();
                site_table.clear();
            }

            void
            record_preserved_nodes(const std::vector<TS_NODE_INT>& node_ids)
            /// Take a list of nodes to record as "ancient samples".
            /// Throws an exception if the nodes are already recorded as such.
            {
                for (auto i : node_ids)
                    {
                        if (static_cast<std::size_t>(i) >= node_table.size())
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

            std::size_t
            push_back_node(double time, std::int32_t pop)
            {
                node_table.push_back(node{ pop, time });
                return node_table.size() - 1;
            }

            template <typename... args>
            std::size_t
            emplace_back_node(args&&... Args)
            {
                node_table.emplace_back(node{ std::forward<args>(Args)... });
                return node_table.size() - 1;
            }

            std::size_t
            push_back_edge(double l, double r, TS_NODE_INT parent,
                           TS_NODE_INT child)
            {
                edge_table.push_back(edge{ l, r, parent, child });
                return edge_table.size();
            }

            template <typename... args>
            std::size_t
            emplace_back_edge(args&&... Args)
            {
                edge_table.emplace_back(edge{ std::forward<args>(Args)... });
                return edge_table.size();
            }

            template <typename... args>
            std::size_t
            push_back_site(args&&... Args)
            {
                site_table.push_back(site{ std::forward<args>(Args)... });
                return site_table.size() - 1;
            }

            template <typename... args>
            std::size_t
            emplace_back_site(args&&... Args)
            {
                site_table.emplace_back(site{ std::forward<args>(Args)... });
                return site_table.size() - 1;
            }

            template <typename... args>
            std::size_t
            push_back_mutation(args&&... Args)
            {
                mutation_table.push_back(
                    mutation_record{ std::forward<args>(Args)... });
                return mutation_table.size() - 1;
            }

            template <typename... args>
            std::size_t
            emplace_back_mutation(args&&... Args)
            {
                mutation_table.emplace_back(
                    mutation_record{ std::forward<args>(Args)... });
                return mutation_table.size() - 1;
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
                        input_left.emplace_back(e.left,
                                                -node_table[e.parent].time,
                                                e.parent, e.child);
                        output_right.emplace_back(e.right,
                                                  node_table[e.parent].time,
                                                  e.parent, e.child);
                    }
                std::sort(input_left.begin(), input_left.end());
                std::sort(output_right.begin(), output_right.end());
            }

            std::size_t
            register_diploid_offspring(
                const std::vector<double>& breakpoints,
                const std::tuple<TS_NODE_INT, TS_NODE_INT>& parents,
                const std::int32_t population, const double time)
            {
                // TODO: carefully document how to index node times.
                auto next_index = emplace_back_node(population, time);
                if (next_index >= std::numeric_limits<TS_NODE_INT>::max())
                    {
                        throw std::invalid_argument("node index too large");
                    }
                split_breakpoints(breakpoints, parents, next_index);
                return next_index;
            }

            void
            rebuild_site_table()
            /// Complete rebuild of the site table.
            {
                auto site_table_copy(site_table);
                site_table.clear();
                for (auto& mr : mutation_table)
                    {
                        auto os = mr.site;
                        record_site_during_rebuild(site_table_copy[mr.site],
                                                   mr);
                        if (site_table_copy[os].position
                            != site_table[mr.site].position)
                            {
                                throw tables_error(
                                    "error rebuilding site table");
                            }
                    }
            }

            void
            sort_mutations_rebuild_site_table()
            /// O(n*log(n)) time plus O(n) additional memory.
            /// This is called by sort_tables, but also should
            /// be called after adding mutations by some manual
            /// means to a table collection.
            {
                sort_mutations();
                rebuild_site_table();
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
            double
            genome_length() const
            {
                return L;
            }
        };

        inline bool
        operator==(const table_collection& a, const table_collection& b)
        {
            return a.genome_length() == b.genome_length()
                   && a.site_table == b.site_table
                   && a.edge_table == b.edge_table
                   && a.node_table == b.node_table
                   && a.mutation_table == b.mutation_table
                   && a.preserved_nodes == b.preserved_nodes;
        }

        inline bool
        operator!=(const table_collection& a, const table_collection& b)
        {
            return !(a == b);
        }

        template <typename mcont_t>
        void
        record_mutations_infinite_sites(
            const TS_NODE_INT u, const mcont_t& mutations,
            const std::vector<std::uint32_t>& new_mutation_keys,
            table_collection& tables)
        /// \version Added in 0.8.0
        {
            std::int8_t ancestral_state = 0, derived_state = 1;
            for (auto& k : new_mutation_keys)
                {
                    auto site = tables.emplace_back_site(mutations[k].pos,
                                                         ancestral_state);
                    if (site >= std::numeric_limits<TS_NODE_INT>::max())
                        {
                            throw std::invalid_argument(
                                "site index out of range");
                        }
                    tables.push_back_mutation(u, k, site, derived_state,
                                              mutations[k].neutral);
                }
        }
    } // namespace ts
} // namespace fwdpp

#endif
