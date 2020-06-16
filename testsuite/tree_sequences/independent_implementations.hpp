/* Independent/haive/slow implementations of tree sequence functions
 * for comparing to efficient versions in the library.
 */
#ifndef FWDPP_TESTUITE_TREE_SEQUENCES_INDEPENDENT_IMPLEMENTATIONS_HPP
#define FWDPP_TESTUITE_TREE_SEQUENCES_INDEPENDENT_IMPLEMENTATIONS_HPP

#include <cstdint>
#include <vector>
#include <fwdpp/ts/node.hpp>
#include <fwdpp/ts/marginal_tree.hpp>

void get_tip(const fwdpp::ts::marginal_tree &m, fwdpp::ts::table_index_t u,
             std::vector<fwdpp::ts::table_index_t> &samples);

std::vector<fwdpp::ts::table_index_t>
naive_get_samples(const fwdpp::ts::marginal_tree &m, fwdpp::ts::table_index_t u);

std::size_t naive_num_samples(const fwdpp::ts::marginal_tree &m,
                              fwdpp::ts::table_index_t u);

// Sample-to-root implementation requiring O(nnodes) extra memory,
// and makes multiple passes through same ancestral nodes.
double
naive_branch_length(const fwdpp::ts::marginal_tree &m,
                    const std::vector<fwdpp::ts::node> &nodes,
                    bool scale_by_length);

#endif
