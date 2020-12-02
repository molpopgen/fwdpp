#include "fwdpp/ts/simplify_tables_output.hpp"
#include <string>
#include <iostream>
#include <numeric>
#include <unordered_set>
#include <boost/program_options.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <gsl/gsl_randist.h>
#include <fwdpp/GSLrng_t.hpp>
#include <fwdpp/ts/std_table_collection.hpp>
#include <fwdpp/ts/table_collection_functions.hpp>
#include <fwdpp/ts/recording.hpp>
#include <fwdpp/ts/recording/diploid_offspring.hpp>
#include <fwdpp/ts/make_simplifier_state.hpp>
#include <fwdpp/ts/simplify_tables.hpp>

#include <tskit.h>

namespace po = boost::program_options;

struct command_line_options
{
    unsigned N;
    double psurvival;
    unsigned nsteps;
    unsigned simplification_interval;
    double rho;
    std::string treefile;
    bool buffer_new_edges;
    bool simplify_from_buffer;
    unsigned seed;

    command_line_options();
};

command_line_options::command_line_options()
    : N{1000}, psurvival{0.}, nsteps{1000},
      simplification_interval{100}, rho{0.}, treefile{"treefile.trees"},
      buffer_new_edges{false}, simplify_from_buffer{false}, seed{42}
{
}

void
validate_cli(const command_line_options& options)
{
    if (options.N == 0)
        {
            throw std::invalid_argument("Population size must be > 0");
        }

    if (options.psurvival < 0. || options.psurvival >= 1.
        || std::isfinite(options.psurvival) == false)
        {
            throw std::invalid_argument("psurvival must be 0.0 <= p < 1.0");
        }

    if (options.rho < 0.0 || std::isfinite(options.rho) == false)
        {
            throw std::invalid_argument("rho must be >= 0.0");
        }

    if (options.treefile.empty())
        {
            throw std::invalid_argument("treefile must not be an empty string");
        }
}

po::options_description
generate_main_options(command_line_options& o)
{
    po::options_description options("Simulation options");
    options.add_options()("help", "Display help");
    options.add_options()("N", po::value<decltype(command_line_options::N)>(&o.N),
                          "Diploid population size. Default = 1000.");

    options.add_options()(
        "psurvival", po::value<decltype(command_line_options::psurvival)>(&o.psurvival),
        "Survival probability. Default = 0.0");
    options.add_options()("nsteps",
                          po::value<decltype(command_line_options::nsteps)>(&o.nsteps),
                          "Number of time steps to evolve. Default = 1000.");
    options.add_options()(
        "simplify",
        po::value<decltype(command_line_options::simplification_interval)>(
            &o.simplification_interval),
        "Time steps between simplifications.  Default = 100.");
    options.add_options()("rho", po::value<decltype(command_line_options::rho)>(&o.rho),
                          "Scaled recombination rate, 4Nr.  Default=0.");
    options.add_options()(
        "treefile", po::value<decltype(command_line_options::treefile)>(&o.treefile),
        "Ouput file name.  Default = treefile.trees");
    options.add_options()(
        "buffer", po::bool_switch(&o.buffer_new_edges),
        "If true, use edge buffering algorithm. If not, sort and simplify. Default = "
        "false");
    options.add_options()("simplify_from_buffer",
                          po::bool_switch(&o.simplify_from_buffer),
                          "If true, and if using the edge buffering algorithm "
                          "(--buffer), then simplify from the buffer");
    options.add_options()("seed",
                          po::value<decltype(command_line_options::seed)>(&o.seed),
                          "Random number seed.  Default = 42.");

    return options;
}

// Simulation types and functions now

struct parent
{
    std::size_t index;
    fwdpp::ts::table_index_t node0, node1;
    parent(std::size_t i, fwdpp::ts::table_index_t n0, fwdpp::ts::table_index_t n1)
        : index(i), node0(n0), node1(n1)
    {
    }
};

struct birth
{
    std::size_t index;
    fwdpp::ts::table_index_t p0node0, p0node1, p1node0, p1node1;
    birth(std::size_t i, const parent& p0, const parent& p1)
        : index(i), p0node0(p0.node0), p0node1(p0.node1), p1node0(p1.node0),
          p1node1(p1.node1)
    {
    }
};

void
deaths_and_parents(const fwdpp::GSLrng_mt& rng, const std::vector<parent>& parents,
                   double psurvival, std::vector<birth>& births)
{
    births.clear();
    for (std::size_t i = 0; i < parents.size(); ++i)
        {
            if (gsl_rng_uniform(rng.get()) > psurvival)
                {
                    std::size_t parent0 = gsl_ran_flat(rng.get(), 0, parents.size());
                    std::size_t parent1 = gsl_ran_flat(rng.get(), 0, parents.size());
                    births.emplace_back(i, parents[parent0], parents[parent1]);
                }
        }
}

void
recombination_breakpoints(const fwdpp::GSLrng_mt& rng, double littler, double maxlen,
                          std::vector<double>& breakpoints)
{
    breakpoints.clear();
    auto nxovers = gsl_ran_poisson(rng.get(), littler);
    for (decltype(nxovers) i = 0; i < nxovers; ++i)
        {
            breakpoints.push_back(gsl_ran_flat(rng.get(), 0., maxlen));
        }
    std::sort(begin(breakpoints), end(breakpoints));
    if (!breakpoints.empty())
        {
            breakpoints.emplace_back(std::numeric_limits<double>::max());
        }
    return;
}

void
generate_births(const fwdpp::GSLrng_mt& rng, const std::vector<birth>& births,
                double littler, std::vector<double>& breakpoints, double birth_time,
                bool buffer_new_edges, fwdpp::ts::edge_buffer& new_edges,
                std::vector<parent>& parents, fwdpp::ts::std_table_collection& tables)
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

template <typename SimplificationState>
void
sort_n_simplify(const std::vector<fwdpp::ts::table_index_t>& samples,
                SimplificationState& state, fwdpp::ts::std_table_collection& tables,
                fwdpp::ts::simplify_tables_output& simplification_output)
{
    auto cmp = fwdpp::ts::get_edge_sort_cmp(tables);
    std::sort(begin(tables.edges), end(tables.edges), cmp);
    fwdpp::ts::simplify_tables(samples, fwdpp::ts::simplification_flags{}, state, tables,
                               simplification_output);
}

template <typename SimplificationState>
void
flush_buffer_n_simplify(
    bool simplify_from_buffer, const std::vector<fwdpp::ts::table_index_t>& samples,
    const std::vector<fwdpp::ts::table_index_t>& alive_at_last_simplification,
    fwdpp::ts::simplify_tables_output& simplification_output,
    fwdpp::ts::edge_buffer& new_edges, SimplificationState& state,
    fwdpp::ts::std_table_collection::edge_table& edge_liftover,
    fwdpp::ts::std_table_collection& tables)
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

using table_collection_ptr
    = std::unique_ptr<tsk_table_collection_t,
                      std::function<void(tsk_table_collection_t*)>>;

void
handle_tskit_return_code(int code)
{
    if (code < 0)
        {
            std::ostringstream o;
            o << tsk_strerror(code);
            throw std::runtime_error(o.str());
        }
}

table_collection_ptr
make_table_collection_ptr(double sequence_length)
{
    table_collection_ptr rv(new tsk_table_collection_t(),
                            [](tsk_table_collection_t* tables) {
                                tsk_table_collection_free(tables);
                                delete tables;
                            });
    int err = tsk_table_collection_init(rv.get(), 0);
    handle_tskit_return_code(err);
    rv->sequence_length = sequence_length;
    if (err != 0)
        {
            throw std::runtime_error("could not initialize tsk_table_collection");
        }
    return rv;
}

void
dump_table_collection_to_tskit(const fwdpp::ts::std_table_collection& tables,
                               const std::string treefile, double forward_time,
                               unsigned N)
{
    auto tskit_tables = make_table_collection_ptr(tables.genome_length());
    unsigned i = 0;
    for (auto& n : tables.nodes)
        {
            // Convert time from forwards to backwards
            // and label the last generation (first 2N nodes)
            // as samples.
            int rv = tsk_node_table_add_row(&tskit_tables->nodes,
                                            (i < 2 * N) ? TSK_NODE_IS_SAMPLE : 0, // flag
                                            -1. * (n.time - forward_time),        // time
                                            TSK_NULL, // population
                                            TSK_NULL, // individual
                                            nullptr,  // metadata
                                            0);       // metadata length
            handle_tskit_return_code(rv);
            ++i;
        }
    for (auto& e : tables.edges)
        {
            auto rv = tsk_edge_table_add_row(&tskit_tables->edges, e.left, e.right,
                                             e.parent, e.child, nullptr, 0);
            handle_tskit_return_code(rv);
        }
    auto rv = tsk_table_collection_build_index(tskit_tables.get(), 0);
    handle_tskit_return_code(rv);
    rv = tsk_table_collection_dump(tskit_tables.get(), treefile.c_str(), 0);
    handle_tskit_return_code(rv);
}

void
simulate(const command_line_options& options)
{
    fwdpp::GSLrng_mt rng(options.seed);
    fwdpp::ts::std_table_collection tables(1.0);
    fwdpp::ts::edge_buffer buffer;
    fwdpp::ts::std_table_collection::edge_table edge_liftover;
    auto simplifier_state = fwdpp::ts::make_simplifier_state(tables);
    std::vector<parent> parents;
    for (unsigned i = 0; i < options.N; ++i)
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
    fwdpp::ts::simplify_tables_output simplification_output;
    bool simplified = false;
    double last_time_simplified = 0; //options.nsteps;
    double littler = options.rho / (4. * static_cast<double>(options.N));
    std::vector<double> breakpoints;
    for (unsigned step = 1; step <= options.nsteps; ++step)
        {
            deaths_and_parents(rng, parents, options.psurvival, births);
            generate_births(rng, births, littler, breakpoints, step,
                            options.buffer_new_edges, buffer, parents, tables);
            if (step % options.simplification_interval == 0.)
                {
                    samples.clear();
                    for (auto& p : parents)
                        {
                            samples.push_back(p.node0);
                            samples.push_back(p.node1);
                        }

                    if (options.buffer_new_edges == false)
                        {
                            sort_n_simplify(samples, simplifier_state, tables,
                                            simplification_output);
                        }
                    else
                        {
                            flush_buffer_n_simplify(
                                options.simplify_from_buffer, samples,
                                alive_at_last_simplification, simplification_output,
                                buffer, simplifier_state, edge_liftover, tables);
                        }
                    simplified = true;
                    last_time_simplified = step;
                    //remap parent nodes
                    for (auto& p : parents)
                        {
                            p.node0 = simplification_output.idmap[p.node0];
                            p.node1 = simplification_output.idmap[p.node1];
                        }
                    if (options.buffer_new_edges == true)
                        {
                            alive_at_last_simplification.clear();
                            for (auto& p : parents)
                                {
                                    alive_at_last_simplification.push_back(p.node0);
                                    alive_at_last_simplification.push_back(p.node1);
                                }
                        }
                }
            else
                {
                    simplified = false;
                }
        }
    if (simplified == false)
        {
            samples.clear();
            for (auto& p : parents)
                {
                    samples.push_back(p.node0);
                    samples.push_back(p.node1);
                }
            if (options.buffer_new_edges == false)
                {
                    sort_n_simplify(samples, simplifier_state, tables,
                                    simplification_output);
                }
            else
                {
                    flush_buffer_n_simplify(options.simplify_from_buffer, samples,
                                            alive_at_last_simplification,
                                            simplification_output, buffer,
                                            simplifier_state, edge_liftover, tables);
                }
        }
    dump_table_collection_to_tskit(tables, options.treefile, options.nsteps, options.N);
}

int
main(int argc, char** argv)
{
    command_line_options options;
    auto cli = generate_main_options(options);
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, cli), vm);
    po::notify(vm);
    validate_cli(options);

    if (vm.count("help"))
        {
            std::cout << cli << '\n';
            std::exit(1);
        }

    simulate(options);
}
