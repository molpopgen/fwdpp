/*! \include wfts.cc
 * Wright-Fisher simulation with tree sequences.
 *
 * See the following paper for background and motivation:
 * Kelleher, Jerome, Kevin Thornton, Jaime Ashander, and Peter Ralph. 2018.
 * “Efficient Pedigree Recording for Fast Population Genetics Simulation.”
 * bioRxiv. https://doi.org/10.1101/248500.
 *
 */

#include <cstdio>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <cassert>
#include <fstream>
#include <string>
#include <fwdpp/ts/table_collection.hpp>
#include <fwdpp/ts/table_simplifier.hpp>
#include <fwdpp/ts/count_mutations.hpp>
#include <fwdpp/ts/recycling.hpp>
#include <fwdpp/ts/tree_visitor.hpp>
#include <fwdpp/ts/mutate_tables.hpp>
#include <fwdpp/ts/generate_data_matrix.hpp>
#include <fwdpp/ts/remove_fixations_from_gametes.hpp>
#include <fwdpp/ts/serialization.hpp>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/sugar/popgenmut.hpp>
#include <fwdpp/sugar/slocuspop.hpp>
#include <fwdpp/util.hpp>
#include <fwdpp/fitness_models.hpp>
#include <fwdpp/extensions/callbacks.hpp>
#include <fwdpp/poisson_xover.hpp>
#include <fwdpp/recbinder.hpp>
#include <boost/program_options.hpp>

#include "simplify_tables.hpp"
#include "evolve_generation_ts.hpp"
#include "calculate_fitnesses.hpp"

namespace po = boost::program_options;
using poptype = fwdpp::slocuspop<fwdpp::popgenmut>;
using GSLrng = fwdpp::GSLrng_t<fwdpp::GSL_RNG_MT19937>;

struct diploid_metadata
{
    std::size_t individual;
    double time, fitness;
    fwdpp::ts::TS_NODE_INT n1, n2;
    diploid_metadata(std::size_t i, double t, double w,
                     fwdpp::ts::TS_NODE_INT a, fwdpp::ts::TS_NODE_INT b)
        : individual(i), time(t), fitness(w), n1(a), n2(b)
    {
    }
};

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

template <typename mcont_t>
void
matrix_runtime_test(const fwdpp::ts::table_collection &tables,
                    const std::vector<fwdpp::ts::TS_NODE_INT> &samples,
                    const mcont_t &mutations,
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

template <typename rng>
std::function<double()>
make_dfe(const fwdpp::uint_t N, const rng &r, const double mean,
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

int
main(int argc, char **argv)
{
    fwdpp::uint_t N, gcint = 100;
    double theta, rho, mean = 0.0, shape = 1, mu,
                       scoeff = std::numeric_limits<double>::quiet_NaN(),
                       dominance = 1.0, scaling = 2.0;
    unsigned seed = 42;
    int ancient_sampling_interval = -1;
    int ancient_sample_size = -1, nsam = 0;
    bool leaf_test = false;
    bool matrix_test = false;
    bool preserve_fixations = false;
    std::string filename, sfsfilename;
    po::options_description options("Simulation options"),
        dfeoptions("Distribution of fitness effects"),
        testing("Testing options");
    // clang-format off
    options.add_options()("help", "Display help")
        ("N", po::value<unsigned>(&N), "Diploid population size")
        ("gc", po::value<unsigned>(&gcint),
        "Simplification interval. Default is 100 generations.")
        ("theta", po::value<double>(&theta), "4Nu")
        ("rho", po::value<double>(&rho), "4Nr")
        ("mu", po::value<double>(&mu), "mutation rate to selected variants")
        ("preserve_fixations",po::bool_switch(&preserve_fixations),"If true, do not count mutations and remove fixations during simulation.  Mutation recycling will proceed via the output of mutation simplification.")
        ("seed", po::value<unsigned>(&seed), "Random number seed. Default is 42")
        ("sampling_interval", po::value<int>(&ancient_sampling_interval), 
         "How often to preserve ancient samples.  Default is -1, which means do not preserve any.")
        ("ansam", po::value<int>(&ancient_sample_size),
         "Sample size (no. diploids) of ancient samples to take at each ancient sampling interval.  Default is -1, and must be reset if sampling_interval is used")
		("sfs", po::value<std::string>(&sfsfilename),"Write the neutral site frequency spectrum of a sample to a file")
		("nsam", po::value<int>(&nsam), "Sample size for the site frequency spectrum.  Default is 0.  Change when using --sfs");
        dfeoptions.add_options()
        ("mean", po::value<double>(&mean), "Mean 2Ns of Gamma distribution of selection coefficients. Default 0.0.")
        ("shape", po::value<double>(&shape), "Shape of Gamma distribution of selection coefficients. Default = 1.")
        ("constant",po::value<double>(&scoeff), "Use a constant DFE with fixed selection coefficient s.\nUsing this over-rides gamma DFE parameters.")
        ("h",po::value<double>(&dominance), "Dominance of selected variants.  Default = 1.0")
        ("scaling",po::value<double>(&scaling), "Fitness model scaling is 1, 1+hs, 1+scaling*s for AA, Aa, and aa genotypes, resp.  Default = 2.0");
        testing.add_options()("leaf_test",po::bool_switch(&leaf_test),"Perform very expensive checking on sample list ranges vs. leaf counts")
        ("matrix_test",po::bool_switch(&matrix_test),"Perform run-time test on generating fwdpp::data_matrix objects and validating the row sums")
		("serialization_test",po::value<std::string>(&filename),"Test round-trip to/from a file");
    // clang-format on
    options.add(dfeoptions);
    options.add(testing);
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, options), vm);
    po::notify(vm);

    if (vm.count("help") || argc == 1)
        {
            std::cout << options << '\n';
            std::exit(1);
        }

    // TODO: need parameter validation
    if (theta < 0. || rho < 0.)
        {
            throw std::invalid_argument("rho and theta must be >= 0.0");
        }
    if (N < 1)
        {
            throw std::invalid_argument("N must be > 0");
        }
    if (gcint < 1)
        {
            throw std::invalid_argument(
                "Simplification (gc) interval must be > 0");
        }
    if (mu < 0)
        {
            throw std::invalid_argument(
                "Mutation rate to selected variants must be >= 0");
        }
    else if (mu > 0)
        {
            if (mean == 0.0 && std::isnan(scoeff))
                {
                    throw std::invalid_argument(
                        "mean selection coefficient cannot be zero");
                }
        }
    if (ancient_sampling_interval > 0 && ancient_sample_size < 1)
        {
            throw std::invalid_argument(
                "ansam must be > 0 when tracking ancient samples");
        }

    GSLrng rng(seed);

    poptype pop(N);
    fwdpp::ts::table_collection tables(2 * pop.diploids.size(), 0, 0, 1.0);
    fwdpp::ts::table_simplifier simplifier(1.0);
    unsigned generation = 1;
    double recrate = rho / static_cast<double>(4 * N);
    const auto recmap
        = fwdpp::recbinder(fwdpp::poisson_xover(recrate, 0., 1.), rng.get());

    auto get_selection_coefficient = make_dfe(N, rng, mean, shape, scoeff);
    const auto generate_mutation_position
        = [&rng]() { return gsl_rng_uniform(rng.get()); };
    const auto generate_h = [dominance]() { return dominance; };
    const auto mmodel = [&pop, &rng, &generation, generate_mutation_position,
                         get_selection_coefficient,
                         generate_h](std::queue<std::size_t> &recbin,
                                     poptype::mcont_t &mutations) {
        return fwdpp::infsites_popgenmut(
            recbin, mutations, rng.get(), pop.mut_lookup, generation,
            // 1.0 signifies 100% of mutations will be selected
            1.0, generate_mutation_position, get_selection_coefficient,
            generate_h);
    };

    // Evolve pop for 20N generations
    fwdpp::ts::TS_NODE_INT first_parental_index = 0,
                           next_index = 2 * pop.diploids.size();
    bool simplified = false;
    std::queue<std::size_t> mutation_recycling_bin;
    std::vector<std::size_t> individual_labels(N);
    std::iota(individual_labels.begin(), individual_labels.end(), 0);
    std::vector<std::size_t> individuals;
    if (ancient_sample_size > 0)
        {
            individuals.resize(ancient_sample_size);
        }
    std::vector<double> fitnesses;
    std::vector<diploid_metadata> ancient_sample_metadata;
    const auto update_offspring = [](std::size_t, std::size_t, std::size_t) {};

    // GOTCHA: we calculate fitnesses and generate our lookup table
    // here, based on our initial/monomorphic population.
    // The reason that we do this is that we will then update
    // these data AFTER each call to the evolution function,
    // so that fitnesses == those of current pop.diploids,
    // meaning that we can record a correct fitness
    // value as ancient sample metadata.  This is a logic
    // issue that is easy to goof.
    auto ff = fwdpp::multiplicative_diploid(scaling);
    auto lookup = calculate_fitnesses(pop, fitnesses, ff);
    for (; generation <= 10 * N; ++generation)
        {
            auto pick1 = [&lookup, &rng]() {
                return gsl_ran_discrete(rng.get(), lookup.get());
            };
            auto pick2 = [&lookup, &rng](const std::size_t /*p1*/) {
                return gsl_ran_discrete(rng.get(), lookup.get());
            };
            evolve_generation(rng, pop, N, mu, pick1, pick2, update_offspring,
                              mmodel, mutation_recycling_bin, recmap,
                              generation, tables, first_parental_index,
                              next_index);
            // Recalculate fitnesses and the lookup table.
            lookup = calculate_fitnesses(pop, fitnesses, ff);
            if (generation % gcint == 0.0)
                {
                    auto rv = simplify_tables(
                        pop, generation, pop.mcounts_from_preserved_nodes, tables,
                        simplifier, tables.num_nodes() - 2 * N, 2 * N,
                        preserve_fixations);
                    if (!preserve_fixations)
                        {
                            mutation_recycling_bin = fwdpp::ts::make_mut_queue(
                                pop.mcounts, pop.mcounts_from_preserved_nodes);
                        }
                    else
                        {
                            mutation_recycling_bin = fwdpp::ts::make_mut_queue(
                                rv.second, pop.mutations.size());
                        }
                    simplified = true;
                    next_index = tables.num_nodes();
                    first_parental_index = 0;

                    // When tracking ancient samples, the node ids of those samples change.
                    // Thus, we need to remap our metadata upon simplification
                    for (auto &md : ancient_sample_metadata)
                        {
                            md.n1 = rv.first[md.n1];
                            md.n2 = rv.first[md.n2];
                            assert(md.n1 != fwdpp::ts::TS_NULL_NODE);
                            assert(md.n2 != fwdpp::ts::TS_NULL_NODE);
                        }
                }
            else
                {
                    simplified = false;
                    first_parental_index = next_index;
                    next_index += 2 * N;

                    // The following (commented-out) block
                    // shows that it is possible to mix mutation
                    // counting strategies in the right situations.
                    // However, the process_gametes call is a quadratic
                    // operation and the overall effect on run time
                    // and peak RAM usage is modest at best in the big
                    // picture...

                    //if (tables.preserved_nodes.empty())
                    //    {
                    //        fwdpp::fwdpp_internal::process_gametes(
                    //            pop.gametes, pop.mutations, pop.mcounts);
                    //        fwdpp::fwdpp_internal::gamete_cleaner(
                    //            pop.gametes, pop.mutations, pop.mcounts, 2 * N,
                    //            std::true_type());
                    //        fwdpp::update_mutations(
                    //            pop.mutations, pop.fixations,
                    //            pop.fixation_times, pop.mut_lookup,
                    //            pop.mcounts, generation, 2 * N);
                    //        mutation_recycling_bin
                    //            = fwdpp::fwdpp_internal::make_mut_queue(
                    //                pop.mcounts);
                    //        auto itr = std::remove_if(
                    //            tables.mutation_table.begin(),
                    //            tables.mutation_table.end(),
                    //            [&pop](const fwdpp::ts::mutation_record &mr) {
                    //                return pop.mcounts[mr.key] == 0;
                    //            });
                    //        tables.mutation_table.erase(
                    //            itr, tables.mutation_table.end());
                    //    }
                }
            if (ancient_sampling_interval > 0
                && generation % ancient_sampling_interval == 0.0
                // This last check forbids us recording the
                // final generation as ancient samples.
                && generation < 10 * N)
                {
                    gsl_ran_choose(
                        rng.get(), individuals.data(), individuals.size(),
                        individual_labels.data(), individual_labels.size(),
                        sizeof(std::size_t));
                    for (auto i : individuals)
                        {
                            auto x = fwdpp::ts::get_parent_ids(
                                first_parental_index, i, 0);
                            assert(x.first >= first_parental_index);
                            assert(x.second >= first_parental_index);
                            assert(x.first < tables.num_nodes());
                            assert(x.second < tables.num_nodes());
                            assert(tables.node_table[x.first].time
                                   == generation);
                            assert(tables.node_table[x.second].time
                                   == generation);
                            assert(std::find(tables.preserved_nodes.begin(),
                                             tables.preserved_nodes.end(),
                                             x.first)
                                   == tables.preserved_nodes.end());
                            assert(std::find(tables.preserved_nodes.begin(),
                                             tables.preserved_nodes.end(),
                                             x.second)
                                   == tables.preserved_nodes.end());
                            tables.preserved_nodes.push_back(x.first);
                            tables.preserved_nodes.push_back(x.second);
                            // Record the metadata for our ancient samples
                            // Here, we record each individual's actual fitness.
                            // If we wanted relative fitness, then
                            // we'd have to copy fitnesses into a temp
                            // and normalize it appropriately.
                            ancient_sample_metadata.emplace_back(
                                i, generation, fitnesses[i], x.first,
                                x.second);
                        }
                }
        }
    if (!simplified)
        {
            auto rv = simplify_tables(pop, generation,
                                      pop.mcounts_from_preserved_nodes, tables,
                                      simplifier, tables.num_nodes() - 2 * N,
                                      2 * N, preserve_fixations);
            if (preserve_fixations)
                {
                    std::vector<std::int32_t> samples(2 * N);
                    std::iota(samples.begin(), samples.end(), 0);
                    fwdpp::ts::count_mutations(tables, pop.mutations, samples,
                                               pop.mcounts,
                                               pop.mcounts_from_preserved_nodes);
                }
            confirm_mutation_counts(pop, tables);
            // When tracking ancient samples, the node ids of those samples change.
            // Thus, we need to remap our metadata upon simplification
            for (auto &md : ancient_sample_metadata)
                {
                    md.n1 = rv.first[md.n1];
                    md.n2 = rv.first[md.n2];
                    assert(md.n1 != fwdpp::ts::TS_NULL_NODE);
                    assert(md.n2 != fwdpp::ts::TS_NULL_NODE);
                }
        }
    // If we have done things correctly, then our
    // ancient sample metadata must match up with
    // what is in our node table.
    for (auto &mr : ancient_sample_metadata)
        {
            if (tables.node_table[mr.n1].time != mr.time
                || tables.node_table[mr.n2].time != mr.time)
                {
                    throw std::runtime_error(
                        "invalid ancient sample metadata");
                }
        }
    assert(tables.input_left.size() == tables.edge_table.size());
    assert(tables.output_right.size() == tables.edge_table.size());
    std::vector<fwdpp::ts::TS_NODE_INT> s(2 * N);
    std::iota(s.begin(), s.end(), 0);
    const auto neutral_variant_maker
        = [&rng, &pop,
           &mutation_recycling_bin](const double left, const double right,
                                    const fwdpp::uint_t generation) {
              return fwdpp::infsites_popgenmut(
                  mutation_recycling_bin, pop.mutations, rng.get(),
                  pop.mut_lookup, generation, 0.0,
                  [left, right, &rng] {
                      return gsl_ran_flat(rng.get(), left, right);
                  },
                  []() { return 0.0; }, []() { return 0.0; });
          };
    auto neutral_muts
        = fwdpp::ts::mutate_tables(rng, neutral_variant_maker, tables, s,
                                   theta / static_cast<double>(4 * N));
    std::sort(tables.mutation_table.begin(), tables.mutation_table.end(),
              [&pop](const fwdpp::ts::mutation_record &a,
                     const fwdpp::ts::mutation_record &b) {
                  return pop.mutations[a.key].pos < pop.mutations[b.key].pos;
              });
    fwdpp::ts::count_mutations(tables, pop.mutations, s, pop.mcounts,
                               pop.mcounts_from_preserved_nodes);
    for (std::size_t i = 0; i < pop.mutations.size(); ++i)
        {
            if (pop.mutations[i].neutral)
                {
                    if (!pop.mcounts[i] && !pop.mcounts_from_preserved_nodes[i])
                        {
                            throw std::runtime_error(
                                "invalid final mutation count");
                        }
                }
        }
    std::cout << neutral_muts << '\n';

    if (leaf_test)
        {
            std::cerr << "Starting sample list validation.  This may take a "
                         "while!\n";
            expensive_leaf_test(tables, s);
            std::cout << "Passed with respect to last generation.\n";
            expensive_leaf_test(tables, tables.preserved_nodes);
            std::cout << "Passed with respect to preserved samples.\n";
        }

    if (matrix_test)
        {
            std::cerr << "Matrix test with respect to last generation...";
            matrix_runtime_test(tables, s, pop.mutations, pop.mcounts);
            std::cerr << "passed\n";
            if (!tables.preserved_nodes.empty())
                {
                    std::cout
                        << "Matrix test with respect to preserved samples...";
                    matrix_runtime_test(tables, tables.preserved_nodes,
                                        pop.mutations,
                                        pop.mcounts_from_preserved_nodes);
                    std::cerr << "passed\n";
                    auto sc = s;
                    sc.insert(sc.end(), tables.preserved_nodes.begin(),
                              tables.preserved_nodes.end());
                    auto mc(pop.mcounts);
                    std::transform(mc.begin(), mc.end(),
                                   pop.mcounts_from_preserved_nodes.begin(),
                                   mc.begin(), std::plus<fwdpp::uint_t>());
                    std::cout << "Matrix test with respect to last generation "
                                 "+ preserved nodes...";
                    matrix_runtime_test(tables, sc, pop.mutations, mc);
                    std::cout << "passed.\n";
                    std::cout << "Matrix test with respect to most recent "
                                 "ancient sampling time point...";
                    sc.clear();
                    std::copy_if(
                        tables.preserved_nodes.begin(),
                        tables.preserved_nodes.end(), std::back_inserter(sc),
                        [&tables](const fwdpp::ts::TS_NODE_INT n) {
                            return tables.node_table[n].time
                                   == tables
                                          .node_table[tables.preserved_nodes
                                                          .back()]
                                          .time;
                        });
                    mc.clear();
                    fwdpp::ts::count_mutations(tables, pop.mutations, sc, mc);
                    matrix_runtime_test(tables, sc, pop.mutations, mc);
                    std::cout << "passed\n";
                }
        }
    if (!filename.empty())
        {
            test_serialization(tables, filename);
        }

    if (!sfsfilename.empty())
        {
            if (!(nsam > 2))
                {
                    throw std::invalid_argument(
                        "sample size for site frequency spectrum must be > 2");
                }
            // Simplify w.r.to 100 samples
            std::vector<fwdpp::ts::TS_NODE_INT> small_sample(nsam);
            gsl_ran_choose(rng.get(), small_sample.data(), small_sample.size(),
                           s.data(), s.size(), sizeof(fwdpp::ts::TS_NODE_INT));
            std::iota(small_sample.begin(), small_sample.end(), 0);
            auto dm = fwdpp::ts::generate_data_matrix(
                tables, small_sample, pop.mutations, true, false);
            auto rs = fwdpp::row_sums(dm);
            std::vector<int> sfs(small_sample.size() - 1);
            for (auto i : rs.first)
                {
                    sfs[i - 1]++;
                }
            std::ofstream sfs_stream(sfsfilename.c_str());
            for (std::size_t i = 0; i < sfs.size(); ++i)
                {
                    sfs_stream << (i + 1) << ' ' << sfs[i] << '\n';
                }
        }
}
