/* Independent/haive/slow implementations of tree sequence functions
 * for comparing to efficient versions in the library.
 */
#ifndef FWDPP_TESTUITE_TREE_SEQUENCES_INDEPENDENT_IMPLEMENTATIONS_HPP
#define FWDPP_TESTUITE_TREE_SEQUENCES_INDEPENDENT_IMPLEMENTATIONS_HPP

#include <cstdint>
#include <vector>
#include <fwdpp/ts/node.hpp>
#include <fwdpp/ts/marginal_tree.hpp>
#include <fwdpp/ts/table_collection.hpp>

void get_tip(const fwdpp::ts::marginal_tree<fwdpp::ts::table_collection::id_type> &m,
             fwdpp::ts::table_collection::id_type u,
             std::vector<fwdpp::ts::table_collection::id_type> &samples);

std::vector<fwdpp::ts::table_collection::id_type> naive_get_samples(
    const fwdpp::ts::marginal_tree<fwdpp::ts::table_collection::id_type> &m,
    fwdpp::ts::table_collection::id_type u);

std::size_t naive_num_samples(
    const fwdpp::ts::marginal_tree<fwdpp::ts::table_collection::id_type> &m,
    fwdpp::ts::table_collection::id_type u);

// Sample-to-root implementation requiring O(nnodes) extra memory,
// and makes multiple passes through same ancestral nodes.
double naive_branch_length(
    const fwdpp::ts::marginal_tree<fwdpp::ts::table_collection::id_type> &m,
    const std::vector<fwdpp::ts::node> &nodes, bool scale_by_length);

#endif
