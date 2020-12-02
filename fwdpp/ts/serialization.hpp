#ifndef FWDPP_TS_SERIALIZATION_HPP
#define FWDPP_TS_SERIALIZATION_HPP

#include <string>
#include <cstdint>
#include <limits>
#include <vector>
#include <stdexcept>
#include <fwdpp/io/scalar_serialization.hpp>
#include "edge.hpp"
#include "node.hpp"
#include "site.hpp"
#include "mutation_record.hpp"
#include "serialization_version.hpp"
#include "table_collection_functions.hpp"

/*! \namespace fwdpp::ts::io
 *  \brief Binary seriazliation of fwdpp::ts::table_collection
 *
 *  The implementation of the binary formats depends on 
 *  template specialiation of function objects for
 *  reading and writing fwdp::ts::node, fwdpp::ts::edge
 *  and fwdpp::ts::mutation_record objects. The templates
 *  are specialized on the format version.  The current
 *  format version is fwdpp::ts::io::TS_TABLES_VERSION.
 *
 *  The function objects themselves are implementation details
 *  that should not be relied upon.  The main functions for library
 *  users in this namespace are:
 *  1. fwdpp::ts::io::serialize_tables
 *  2. fwdpp::ts::io::deserialize_tables
 *
 *  \version 0.7.0 Added to library
 *  \version 0.9.0 Added typename TableCollectionType
 */

namespace fwdpp
{
    namespace ts
    {
        namespace io
        {
            template <std::size_t version> struct serialize_node
            {
                template <typename ostreamtype>
                inline void
                operator()(ostreamtype& o, const node& n) const
                {
                    throw std::runtime_error("invalid serialization version");
                }
            };

            template <> struct serialize_node<TS_TABLES_VERSION>
            {
                template <typename ostreamtype>
                inline void
                operator()(ostreamtype& o, const node& n) const
                {
                    fwdpp::io::scalar_writer sw;
                    sw(o, &n.deme);
                    sw(o, &n.time);
                }
            };

            template <std::size_t version> struct deserialize_node
            {
                template <typename istreamtype>
                inline node
                operator()(istreamtype&) const
                {
                    throw std::runtime_error("invalid serialization version");
                }
            };

            template <> struct deserialize_node<TS_TABLES_VERSION>
            {
                template <typename istreamtype>
                inline node
                operator()(istreamtype& o) const
                {
                    fwdpp::io::scalar_reader sr;
                    table_index_t deme;
                    double time;
                    sr(o, &deme);
                    sr(o, &time);
                    return node{deme, time};
                }
            };

            template <std::size_t version> struct serialize_edge
            {
                template <typename ostreamtype>
                inline void
                operator()(ostreamtype& o, const edge& e) const
                {
                    throw std::runtime_error("invalid serialization version");
                }
            };

            template <> struct serialize_edge<TS_TABLES_VERSION>
            {
                template <typename ostreamtype>
                inline void
                operator()(ostreamtype& o, const edge& e) const
                {
                    fwdpp::io::scalar_writer sw;
                    sw(o, &e.left);
                    sw(o, &e.right);
                    sw(o, &e.parent);
                    sw(o, &e.child);
                }
            };

            template <std::size_t version> struct deserialize_edge
            {
                template <typename istreamtype>
                inline edge
                operator()(istreamtype&) const
                {
                    throw std::runtime_error("invalid serialization version");
                }
            };

            template <> struct deserialize_edge<TS_TABLES_VERSION>
            {
                template <typename istreamtype>
                inline edge
                operator()(istreamtype& i) const
                {
                    fwdpp::io::scalar_reader sr;
                    double left, right;
                    table_index_t parent, child;
                    sr(i, &left);
                    sr(i, &right);
                    sr(i, &parent);
                    sr(i, &child);
                    return edge{left, right, parent, child};
                }
            };

            template <std::size_t version> struct serialize_mutation_record
            {
                template <typename ostreamtype>
                inline void
                operator()(ostreamtype& o, const mutation_record& mr) const
                {
                    throw std::runtime_error("invalid serialization version");
                }
            };

            template <> struct serialize_mutation_record<2>
            /// \version 0.8.0 Changes from TS_TABLES_VERSION to 2
            {
                template <typename ostreamtype>
                inline void
                operator()(ostreamtype& o, const mutation_record& mr) const
                {
                    fwdpp::io::scalar_writer sw;
                    sw(o, &mr.node);
                    sw(o, &mr.key);
                }
            };

            template <std::size_t version> struct deserialize_mutation_record
            {
                template <typename istreamtype>
                inline mutation_record
                operator()(istreamtype&) const
                {
                    throw std::runtime_error("invalid serialization version");
                }
            };

            template <> struct deserialize_mutation_record<2>
            /// \version 0.8.0 Changed from TS_TABLES_VERSION to 2
            /// \note The fields site, derived_state, and neutral
            ///       are filled in with "junk" values.
            {
                template <typename istreamtype>
                inline mutation_record
                operator()(istreamtype& i) const
                {
                    fwdpp::io::scalar_reader sr;
                    table_index_t node;
                    decltype(mutation_record::key) key;
                    sr(i, &node);
                    sr(i, &key);
                    return mutation_record{
                        node, key, std::numeric_limits<std::size_t>::max(),
                        std::numeric_limits<std::int8_t>::max(), true};
                }
            };

            namespace backwards_compat
            {
                struct mutation_record_V2
                /// This was the mutation_record format
                /// prior to TS_TABLES_VERSION 3.
                {
                    std::int32_t node;
                    std::size_t key;
                };
            } // namespace backwards_compat

            template <typename TableCollectionType, typename ostreamtype>
            void
            serialize_tables(ostreamtype& o, const TableCollectionType& tables)
            /*! \brief Write a fwdpp::ts::table_collection to a binary format.
			 *
			 *  \param o A model of std::ostream
			 *  \param tables A fwdpp::ts::table_collection
			 *
			 *  The result is that the data in \a tables are written
			 *  in a machine-readable format to \a o.
			 *
			 *  This function always outputs to the format version
			 *  specified by fwdpp::ts::io::TS_TABLES_VERSION.
             *
             *  \version 0.7.4 Tables are written in a single call
			 */
            {
                o << "fwdppts";
                fwdpp::io::scalar_writer sw;
                sw(o, &TS_TABLES_VERSION);
                auto L = tables.genome_length();
                sw(o, &L);
                sw(o, &tables.edge_offset);
                std::size_t num_edges = tables.edges.size(),
                            num_nodes = tables.num_nodes(),
                            num_mutations = tables.mutations.size(),
                            num_sites = tables.sites.size();
                sw(o, &num_edges);
                sw(o, &num_nodes);
                sw(o, &num_mutations);
                sw(o, &num_sites);
                if (!tables.edges.empty())
                    {
                        o.write(reinterpret_cast<const char*>(tables.edges.data()),
                                tables.edges.size() * sizeof(edge));
                    }
                if (!tables.nodes.empty())
                    {
                        o.write(reinterpret_cast<const char*>(tables.nodes.data()),
                                tables.nodes.size() * sizeof(node));
                    }
                if (!tables.mutations.empty())
                    {
                        o.write(reinterpret_cast<const char*>(tables.mutations.data()),
                                tables.mutations.size() * sizeof(mutation_record));
                    }
                if (!tables.sites.empty())
                    {
                        o.write(reinterpret_cast<const char*>(tables.sites.data()),
                                tables.sites.size() * sizeof(site));
                    }
            }

            template <typename TableCollectionType> struct deserialize_tables
            {
                template <typename istreamtype>
                inline TableCollectionType
                operator()(istreamtype& i)
                /*! \brief Read a fwdpp::ts::table_collection in from a binary stream
			 *
			 *  \param i A model of std::istream
			 *
			 *  \return fwdpp::ts::table_collection
			 *
			 *  \note The return value has its index vectors populated.
			 *  See fwdpp::ts::table_collection::build_indexes
             *
             *  \version 0.7.4 Tables are read in a single call
			 */
                {
                    //Reading data back in has to manage versions
                    char fwdppts[7];
                    i.read(fwdppts, 7);
                    if (std::string(fwdppts, fwdppts + 7) != "fwdppts")
                        {
                            throw std::runtime_error(
                                "input stream is not at the beginning of "
                                "table_collection");
                        }
                    std::uint32_t format;
                    fwdpp::io::scalar_reader sr;
                    sr(i, &format);
                    double L;
                    sr(i, &L);
                    TableCollectionType tables(L);
                    sr(i, &tables.edge_offset);
                    std::size_t num_edges, num_nodes, num_mutations, num_sites;
                    sr(i, &num_edges);
                    sr(i, &num_nodes);
                    sr(i, &num_mutations);
                    if (format == TS_TABLES_VERSION)
                        {
                            sr(i, &num_sites);
                        }
                    if (format == TS_TABLES_VERSION || format == 2)
                        {
                            tables.edges.resize(num_edges);
                            i.read(reinterpret_cast<char*>(tables.edges.data()),
                                   num_edges * sizeof(edge));
                            tables.nodes.resize(num_nodes);
                            i.read(reinterpret_cast<char*>(tables.nodes.data()),
                                   num_nodes * sizeof(node));

                            if (format == TS_TABLES_VERSION)
                                {
                                    tables.mutations.resize(num_mutations);
                                    i.read(
                                        reinterpret_cast<char*>(tables.mutations.data()),
                                        num_mutations * sizeof(mutation_record));
                                    tables.sites.resize(num_sites);
                                    i.read(reinterpret_cast<char*>(tables.sites.data()),
                                           num_sites * sizeof(site));
                                }
                            else
                                {
                                    std::vector<backwards_compat::mutation_record_V2>
                                        temp(num_mutations);
                                    i.read(
                                        reinterpret_cast<char*>(temp.data()),
                                        num_mutations
                                            * sizeof(
                                                backwards_compat::mutation_record_V2));
                                    for (auto t : temp)
                                        {
                                            tables.mutations.emplace_back(
                                                mutation_record{t.node, t.key,
                                                                std::numeric_limits<
                                                                    std::size_t>::max(),
                                                                std::numeric_limits<
                                                                    std::int8_t>::min(),
                                                                true});
                                        }
                                }
                        }
                    else if (format == 1)
                        {
                            deserialize_edge<TS_TABLES_VERSION> edge_reader;
                            deserialize_node<TS_TABLES_VERSION> node_reader;
                            // Format versions 1 and 2 have the same mutation_record type
                            deserialize_mutation_record<2> mutation_record_reader;
                            tables.edges.reserve(num_edges);
                            for (std::size_t j = 0; j < num_edges; ++j)
                                {
                                    tables.edges.emplace_back(edge_reader(i));
                                }
                            tables.nodes.reserve(num_nodes);
                            for (std::size_t j = 0; j < num_nodes; ++j)
                                {
                                    tables.nodes.emplace_back(node_reader(i));
                                }
                            tables.mutations.reserve(num_mutations);
                            for (std::size_t j = 0; j < num_mutations; ++j)
                                {
                                    tables.mutations.emplace_back(
                                        mutation_record_reader(i));
                                }
                        }
                    else
                        {
                            throw std::runtime_error(
                                "invalid serialization version detected");
                        }
                    if (TS_TABLES_VERSION < 4)
                        {
                            std::size_t num_preserved_samples;
                            sr(i, &num_preserved_samples);
                            if (num_preserved_samples)
                                {
                                    std::vector<table_index_t> preserved_nodes;
                                    preserved_nodes.resize(num_preserved_samples);
                                    sr(i, preserved_nodes.data(), num_preserved_samples);
                                }
                        }
                    tables.build_indexes();
                    return tables;
                }
            };

            template <typename TableCollectionType, typename MutationContainerType>
            inline void
            fix_mutation_table_repopulate_site_table(
                TableCollectionType& tables, const MutationContainerType& mutations)
            /// \brief Helper function when reading back from old file formats
            /// \version 0.8.0 Added to library
            ///
            /// When calling deserialize_tables to read in table_collection
            /// objects generated with fwdpp versions < 0.8.0, the
            /// mutation_record object differs from the current version
            /// and no site_table is present.  This function fixes that.
            {
                tables.sites.clear();
                for (auto& mr : tables.mutations)
                    {
                        if (mr.key >= mutations.size())
                            {
                                throw std::runtime_error("mutation key out of range");
                            }
                        mr.neutral = mutations[mr.key].neutral;
                        mr.derived_state = default_derived_state;
                        mr.site = tables.emplace_back_site(mutations[mr.key].pos,
                                                           default_ancestral_state);
                    }

                // This is almost certainly not necessary,
                // but we may as well be super-extra safe.
                // Calling this will fix any issues due
                // to recording the same position more than
                // once in the site table.
                sort_mutation_table_and_rebuild_site_table(tables);
            }
        } // namespace io
    }     // namespace ts
} // namespace fwdpp

#endif
