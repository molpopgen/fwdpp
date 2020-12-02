#ifndef FWDPP_TESTSUITE_TREE_SEQUENCES_WFEVOLVE_TABLE_COLLECTION_HPP
#define FWDPP_TESTSUITE_TREE_SEQUENCES_WFEVOLVE_TABLE_COLLECTION_HPP

#include <vector>
#include <cstdint>
#include <numeric>
#include <algorithm>
#include <fwdpp/GSLrng_t.hpp>
#include <fwdpp/ts/definitions.hpp>
#include <fwdpp/ts/make_simplifier_state.hpp>
#include <fwdpp/ts/simplify_tables.hpp>
#include <fwdpp/ts/table_collection_functions.hpp>
#include <fwdpp/ts/recording/edge_buffer.hpp>
#include <fwdpp/ts/recording/diploid_offspring.hpp>

#include <gsl/gsl_randist.h>

struct parent
{
    std::size_t index;
    fwdpp::ts::table_index_t nodes[2];
    parent(std::size_t i, fwdpp::ts::table_index_t n0, fwdpp::ts::table_index_t n1)
        : index(i), nodes{n0, n1}
    {
    }
};

struct birth
{
    std::size_t index;
    fwdpp::ts::table_index_t p0node0, p0node1, p1node0, p1node1;
    birth(std::size_t i, const parent& p0, const parent& p1)
        : index(i), p0node0(p0.nodes[0]), p0node1(p0.nodes[1]), p1node0(p1.nodes[0]),
          p1node1(p1.nodes[1])
    {
    }
};

struct wf_simulation_params
{
    unsigned seed;
    unsigned N;
    unsigned nsteps;
    double psurvival;
    double rho;
    unsigned simplification_interval;
};

struct empty_policies
{
};

struct wf_simulation_results
// The return value of the simulation
{
    std::vector<parent> alive_individuals;
    std::vector<fwdpp::ts::table_index_t> preserved_nodes;

    template <typename A, typename P>
    wf_simulation_results(A&& a, P&& p)
        : alive_individuals{std::forward<A>(a)}, preserved_nodes{std::forward<P>(p)}
    {
    }
};

void deaths_and_parents(const fwdpp::GSLrng_mt& rng, const std::vector<parent>& parents,
                        double psurvival, std::vector<birth>& births);

void recombination_breakpoints(const fwdpp::GSLrng_mt& rng, double littler,
                               double maxlen, std::vector<double>& breakpoints);

template <typename TableCollectionType>
inline void
generate_births(const fwdpp::GSLrng_mt& rng, const std::vector<birth>& births,
                double littler, std::vector<double>& breakpoints, double birth_time,
                bool buffer_new_edges, fwdpp::ts::edge_buffer& new_edges,
                std::vector<parent>& parents, TableCollectionType& tables)
{
    std::size_t new_node_0, new_node_1;
    for (auto& b : births)
        {
            auto p0n0 = b.p0node0;
            auto p0n1 = b.p0node1;
            if (gsl_rng_uniform(rng.get()) < 0.5)
                {
                    std::swap(p0n0, p0n1);
                }
            auto p1n0 = b.p1node0;
            auto p1n1 = b.p1node1;
            if (gsl_rng_uniform(rng.get()) < 0.5)
                {
                    std::swap(p1n0, p1n1);
                }
            if (buffer_new_edges == false)
                {
                    recombination_breakpoints(rng, littler, tables.genome_length(),
                                              breakpoints);
                    new_node_0 = fwdpp::ts::record_diploid_offspring(
                        breakpoints, std::tie(p0n0, p0n1), 0, birth_time, tables);
                    recombination_breakpoints(rng, littler, tables.genome_length(),
                                              breakpoints);
                    new_node_1 = fwdpp::ts::record_diploid_offspring(
                        breakpoints, std::tie(p1n0, p1n1), 0, birth_time, tables);
                }
            else
                {
                    recombination_breakpoints(rng, littler, tables.genome_length(),
                                              breakpoints);
                    new_node_0 = fwdpp::ts::record_diploid_offspring(
                        breakpoints, std::tie(p0n0, p0n1), 0, birth_time, tables,
                        new_edges);
                    double ptime = tables.nodes[p0n0].time;
                    double ctime = tables.nodes[new_node_0].time;
                    if (ctime <= ptime)
                        {
                            throw std::runtime_error("bad parent/child time");
                        }
                    recombination_breakpoints(rng, littler, tables.genome_length(),
                                              breakpoints);
                    new_node_1 = fwdpp::ts::record_diploid_offspring(
                        breakpoints, std::tie(p1n0, p1n1), 0, birth_time, tables,
                        new_edges);
                    ptime = tables.nodes[p1n0].time;
                    ctime = tables.nodes[new_node_1].time;
                    if (ctime <= ptime)
                        {
                            throw std::runtime_error("bad parent/child time");
                        }
                }
            parents[b.index] = parent(b.index, new_node_0, new_node_1);
        }
}

template <typename TableCollectionType, typename SimplificationState>
void
sort_n_simplify(const std::vector<fwdpp::ts::table_index_t>& samples,
                SimplificationState& state, TableCollectionType& tables,
                fwdpp::ts::simplify_tables_output& simplification_output)
{
    auto cmp = fwdpp::ts::get_edge_sort_cmp(tables);
    std::sort(begin(tables.edges), end(tables.edges), cmp);
    fwdpp::ts::simplify_tables(samples, fwdpp::ts::simplification_flags{}, state, tables,
                               simplification_output);
}

template <typename TableCollectionType, typename SimplificationState>
void
flush_buffer_n_simplify(
    bool simplify_from_buffer, const std::vector<fwdpp::ts::table_index_t>& samples,
    const std::vector<fwdpp::ts::table_index_t>& alive_at_last_simplification,
    fwdpp::ts::simplify_tables_output& simplification_output,
    fwdpp::ts::edge_buffer& new_edges, SimplificationState& state,
    typename TableCollectionType::edge_table& edge_liftover, TableCollectionType& tables)
{
    double max_time = -1; //-1;//std::numeric_limits<double>::max();
    for (auto a : alive_at_last_simplification)
        {
            max_time = std::max(max_time, tables.nodes[a].time);
        }

    if (simplify_from_buffer == false)
        {
            stitch_together_edges(alive_at_last_simplification, max_time, new_edges,
                                  edge_liftover, tables);
            fwdpp::ts::simplify_tables(samples, fwdpp::ts::simplification_flags{}, state,
                                       tables, simplification_output);
            // stitch_together_edges doesn't know about what happens during
            // simplify, so we need to manually reset the buffer's head/tail
            // sizes
            new_edges.reset(tables.num_nodes());
        }
    else
        {
            // Simplifying right from the buffer automatically
            // resets new_edges
            fwdpp::ts::simplify_tables(samples, alive_at_last_simplification,
                                       fwdpp::ts::simplification_flags{}, state, tables,
                                       new_edges, simplification_output);
        }
}

// NOTE: the SimulationPolicies type is a placeholder
// for functions to do things like add mutations, remember
// samples, etc., etc..  Future updates will probably have
// to metaprogram to make sure that empty structs can be used.
template <typename TableCollectionType, typename SimulationPolicies>
wf_simulation_results
wfevolve_table_collection(unsigned seed, unsigned N, unsigned nsteps, double psurvival,
                          double rho, unsigned simplification_interval,
                          bool buffer_new_edges, bool simplify_from_buffer,
                          bool never_simplify, SimulationPolicies /*policies*/,
                          TableCollectionType& tables)
{
    fwdpp::GSLrng_mt rng(seed);
    fwdpp::ts::edge_buffer buffer;
    fwdpp::ts::simplify_tables_output simplification_output;
    typename TableCollectionType::edge_table edge_liftover;
    auto simplifier_state = fwdpp::ts::make_simplifier_state(tables);
    std::vector<parent> parents;
    for (unsigned i = 0; i < N; ++i)
        {
            auto id0 = tables.emplace_back_node(0, 0.);
            auto id1 = tables.emplace_back_node(0, 0.);
            parents.emplace_back(i, id0, id1);
        }

    // The next bits are all for buffering
    std::vector<fwdpp::ts::table_index_t> alive_at_last_simplification(
        tables.num_nodes());
    std::iota(begin(alive_at_last_simplification), end(alive_at_last_simplification), 0);

    std::vector<birth> births;
    std::vector<fwdpp::ts::table_index_t> samples;
    bool simplified = false;
    double littler = rho / (4. * static_cast<double>(N));

    std::vector<double> breakpoints;
    for (unsigned step = 1; step <= nsteps; ++step)
        {
            deaths_and_parents(rng, parents, psurvival, births);
            generate_births(rng, births, littler, breakpoints, step, buffer_new_edges,
                            buffer, parents, tables);
            if (never_simplify == false && step % simplification_interval == 0.)
                {
                    samples.clear();
                    for (auto& p : parents)
                        {
                            samples.push_back(p.nodes[0]);
                            samples.push_back(p.nodes[1]);
                        }

                    if (buffer_new_edges == false)
                        {
                            sort_n_simplify(samples, simplifier_state, tables,
                                            simplification_output);
                        }
                    else
                        {
                            flush_buffer_n_simplify(
                                simplify_from_buffer, samples,
                                alive_at_last_simplification, simplification_output,
                                buffer, simplifier_state, edge_liftover, tables);
                        }
                    simplified = true;
                    //remap parent nodes
                    for (auto& p : parents)
                        {
                            p.nodes[0] = simplification_output.idmap[p.nodes[0]];
                            p.nodes[1] = simplification_output.idmap[p.nodes[1]];
                        }
                    if (buffer_new_edges == true)
                        {
                            alive_at_last_simplification.clear();
                            for (auto& p : parents)
                                {
                                    alive_at_last_simplification.push_back(p.nodes[0]);
                                    alive_at_last_simplification.push_back(p.nodes[1]);
                                }
                        }
                }
            else
                {
                    simplified = false;
                }
        }
    if (never_simplify == false && simplified == false)
        {
            samples.clear();
            for (auto& p : parents)
                {
                    samples.push_back(p.nodes[0]);
                    samples.push_back(p.nodes[1]);
                }
            if (buffer_new_edges == false)
                {
                    sort_n_simplify(samples, simplifier_state, tables,
                                    simplification_output);
                }
            else
                {
                    flush_buffer_n_simplify(simplify_from_buffer, samples,
                                            alive_at_last_simplification,
                                            simplification_output, buffer,
                                            simplifier_state, edge_liftover, tables);
                }
        }
    if (never_simplify == true && buffer_new_edges == true)
        {
            // Then we've kinda wasted some effort and we need to
            // copy stuff to the edge table. ;)
            double max_time = -1; //-1;//std::numeric_limits<double>::max();
            for (auto a : alive_at_last_simplification)
                {
                    max_time = std::max(max_time, tables.nodes[a].time);
                }
            stitch_together_edges(alive_at_last_simplification, max_time, buffer,
                                  edge_liftover, tables);
        }

    return wf_simulation_results{std::move(parents),
                                 std::vector<fwdpp::ts::table_index_t>{}};
}

#endif
