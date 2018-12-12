#include <cmath>
#include <stdexcept>
#include <fstream>
#include <fwdpp/ts/generate_data_matrix.hpp>
#include <fwdpp/ts/serialization.hpp>
#include <fwdpp/extensions/callbacks.hpp>
#include "tree_sequence_examples_common.hpp"

std::function<double()>
make_dfe(const fwdpp::uint_t N, const fwdpp::GSLrng_mt &r, const double mean,
         const double shape, const double scoeff)
{
    if (std::isfinite(scoeff))
        {
            return [scoeff]() { return scoeff; };
        }
    fwdpp::extensions::gamma dfe(mean, shape);
    return
        [&r, dfe, N]() { return dfe(r.get()) / static_cast<double>(2 * N); };
}

void
matrix_runtime_test(const fwdpp::ts::table_collection &tables,
                    const std::vector<fwdpp::ts::TS_NODE_INT> &samples,
                    const std::vector<fwdpp::popgenmut> &mutations,
                    const std::vector<fwdpp::uint_t> &mcounts)
{
    auto dm = fwdpp::ts::generate_data_matrix(tables, samples, mutations, true,
                                              true);
    auto rs = fwdpp::row_sums(dm);
    for (std::size_t i = 0; i < rs.first.size(); ++i)
        {
            if (rs.first[i] != mcounts[dm.neutral_keys[i]])
                {
                    throw std::runtime_error("bad neutral mutation count");
                }
        }
    for (std::size_t i = 0; i < rs.second.size(); ++i)
        {
            if (rs.second[i] != mcounts[dm.selected_keys[i]])
                {
                    throw std::runtime_error("bad selected mutation count");
                }
        }
}

void
expensive_leaf_test(const fwdpp::ts::table_collection &tables,
                    const std::vector<fwdpp::ts::TS_NODE_INT> &sample_list)
{
    fwdpp::ts::tree_visitor mti(tables, sample_list);
    while (mti(std::true_type(), std::true_type()))
        {
            auto &tree = mti.tree();
            for (auto i : sample_list)
                {
                    auto p = i;
                    while (p != -1)
                        {
                            auto l = tree.left_sample[p];
                            auto ogl = l;
                            if (l != -1)
                                {
                                    auto r = tree.right_sample[p];
                                    int ns = 0;
                                    while (true)
                                        {
                                            ++ns;
                                            if (l == r)
                                                {
                                                    break;
                                                }
                                            l = tree.next_sample[l];
                                            if (l == ogl)
                                                {
                                                    throw std::runtime_error(
                                                        "loopback error");
                                                }
                                        }
                                    if (ns != tree.leaf_counts[p])
                                        {
                                            throw std::runtime_error(
                                                "bad sample interval");
                                        }
                                }
                            p = tree.parents[p];
                        }
                }
        }
}

void
test_serialization(const fwdpp::ts::table_collection &tables,
                   const std::string &filename)
{
    std::ofstream o(filename.c_str());
    fwdpp::ts::io::serialize_tables(o, tables);
    o.close();
    std::ifstream i(filename.c_str());
    auto tables2 = fwdpp::ts::io::deserialize_tables(i);

    if (tables.genome_length() != tables2.genome_length())
        {
            throw std::runtime_error("genome_length does not match");
        }
    if (tables.edge_offset != tables2.edge_offset)
        {
            throw std::runtime_error("edge_offset does not match");
        }

    if (tables.edge_table != tables2.edge_table)
        {
            throw std::runtime_error("edge tables do not match");
        }

    if (tables.node_table != tables2.node_table)
        {
            throw std::runtime_error("node tables do not match");
        }

    if (tables.mutation_table != tables2.mutation_table)
        {
            throw std::runtime_error("mutation tables do not match");
        }
    if (tables != tables2)
        {
            throw std::runtime_error("tables failed equality check");
        }
}
