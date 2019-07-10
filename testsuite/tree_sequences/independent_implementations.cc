#include <algorithm>
#include <fwdpp/ts/node.hpp>
#include <fwdpp/ts/marginal_tree.hpp>
#include <fwdpp/ts/marginal_tree_functions/children.hpp>

void
get_tip(const fwdpp::ts::marginal_tree &m, fwdpp::ts::TS_NODE_INT u,
        std::vector<fwdpp::ts::TS_NODE_INT> &samples)
{
    if (fwdpp::ts::num_children(m, u) > 0)
        {
            fwdpp::ts::process_children(
                m, u, true, [&m, &samples](fwdpp::ts::TS_NODE_INT x) {
                    get_tip(m, x, samples);
                });
            return;
        }
    samples.push_back(u);
}

std::vector<fwdpp::ts::TS_NODE_INT>
naive_get_samples(const fwdpp::ts::marginal_tree &m, fwdpp::ts::TS_NODE_INT u)
{
    if (!m.advancing_sample_list())
        {
            throw std::invalid_argument("sample lists are not being updated");
        }
    std::vector<fwdpp::ts::TS_NODE_INT> temp,
        samples_list(m.samples_list_begin(), m.samples_list_end()),
        intersection;
    std::sort(begin(samples_list), end(samples_list));
    get_tip(m, u, temp);
    std::sort(begin(temp), end(temp));
    std::set_intersection(begin(temp), end(temp), begin(samples_list),
                          end(samples_list), std::back_inserter(intersection));
    return intersection;
}

std::size_t
naive_num_samples(const fwdpp::ts::marginal_tree &m, fwdpp::ts::TS_NODE_INT u)
{
    auto s = naive_get_samples(m, u);
    return s.size();
}

double
naive_branch_length(const fwdpp::ts::marginal_tree &m,
                    const std::vector<fwdpp::ts::node> &nodes,
                    bool scale_by_length)
// Sample-to-root implementation requiring O(nnodes) extra memory,
// and makes multiple passes through same ancestral nodes.
{
    double ttime = 0;
    std::vector<std::int8_t> processed(nodes.size(), 0);
    for (auto s = m.samples_list_begin(); s != m.samples_list_end(); ++s)
        {
            auto u = *s;
            while (u != fwdpp::ts::TS_NULL_NODE)
                {
                    auto p = m.parents[u];
                    if (!processed[u] && p != fwdpp::ts::TS_NULL_NODE)
                        {
                            ttime += (nodes[u].time - nodes[p].time);
                        }
                    processed[u] = 1;
                    u = p;
                }
        }
    if (scale_by_length)
        {
            ttime *= (m.right - m.left);
        }
    return ttime;
}
