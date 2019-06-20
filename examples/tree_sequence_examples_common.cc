#include <cmath>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <gsl/gsl_randist.h>
#include <fwdpp/ts/count_mutations.hpp>
#include <fwdpp/ts/generate_data_matrix.hpp>
#include <fwdpp/ts/serialization.hpp>
#include <fwdpp/ts/mutate_tables.hpp>
#include <fwdpp/extensions/callbacks.hpp>
#include "tree_sequence_examples_common.hpp"

namespace po = boost::program_options;

template <typename poptype>
int
apply_neutral_mutations_details(
    const options &o, const fwdpp::GSLrng_mt &rng,
    fwdpp::ts::table_collection &tables, poptype &pop,
    fwdpp::flagged_mutation_queue &mutation_recycling_bin)
{
    std::vector<fwdpp::ts::TS_NODE_INT> s(2 * o.N);

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
                                   o.theta / static_cast<double>(4 * o.N));
    std::sort(tables.mutation_table.begin(), tables.mutation_table.end(),
              [&pop](const fwdpp::ts::mutation_record &a,
                     const fwdpp::ts::mutation_record &b) {
                  return pop.mutations[a.key].pos < pop.mutations[b.key].pos;
              });
    return neutral_muts;
}

int
apply_neutral_mutations(const options &o, const fwdpp::GSLrng_mt &rng,
                        fwdpp::ts::table_collection &tables,
                        single_locus_poptype &pop,
                        fwdpp::flagged_mutation_queue &mutation_recycling_bin)
{
    return apply_neutral_mutations_details(o, rng, tables, pop,
                                           mutation_recycling_bin);
}

options::options()
    : N{}, gcint(100), theta(), rho(), mean(0.), shape(1.), mu(),
      scoeff(std::numeric_limits<double>::max()), dominance(1.), scaling(2.),
      seed(42), ancient_sampling_interval(-1), ancient_sample_size(-1),
      nsam(0), leaf_test(false), matrix_test(false),
      preserve_fixations(false), filename(), sfsfilename()
{
}

po::options_description
generate_main_options(options &o)
{
    po::options_description options("Simulation options");
    // clang-format off
    options.add_options()("help", "Display help")
        ("N", po::value<unsigned>(&o.N), "Diploid population size")
        ("gc", po::value<unsigned>(&o.gcint),
        "Simplification interval. Default is 100 generations.")
        ("theta", po::value<double>(&o.theta), "4Nu")
        ("rho", po::value<double>(&o.rho), "4Nr")
        ("mu", po::value<double>(&o.mu), "mutation rate to selected variants")
        ("preserve_fixations",po::bool_switch(&o.preserve_fixations),"If true, do not count mutations and remove fixations during simulation.  Mutation recycling will proceed via the output of mutation simplification.")
        ("seed", po::value<unsigned>(&o.seed), "Random number seed. Default is 42")
        ("sampling_interval", po::value<int>(&o.ancient_sampling_interval), 
         "How often to preserve ancient samples.  Default is -1, which means do not preserve any.")
        ("ansam", po::value<int>(&o.ancient_sample_size),
         "Sample size (no. diploids) of ancient samples to take at each ancient sampling interval.  Default is -1, and must be reset if sampling_interval is used")
		("sfs", po::value<std::string>(&o.sfsfilename),"Write the neutral site frequency spectrum of a sample to a file")
		("nsam", po::value<int>(&o.nsam), "Sample size for the site frequency spectrum.  Default is 0.  Change when using --sfs");
    // clang-format on
    return options;
}

po::options_description
generate_dfe_options(options &o)
{
    po::options_description dfeoptions("Distribution of fitness effects");
    // clang-format off
    dfeoptions.add_options()
        ("mean", po::value<double>(&o.mean), "Mean 2Ns of Gamma distribution of selection coefficients. Default 0.0.")
        ("shape", po::value<double>(&o.shape), "Shape of Gamma distribution of selection coefficients. Default = 1.")
        ("constant",po::value<double>(&o.scoeff), "Use a constant DFE with fixed selection coefficient s.\nUsing this over-rides gamma DFE parameters.")
        ("h",po::value<double>(&o.dominance), "Dominance of selected variants.  Default = 1.0")
        ("scaling",po::value<double>(&o.scaling), "Fitness model scaling is 1, 1+hs, 1+scaling*s for AA, Aa, and aa genotypes, resp.  Default = 2.0");
    // clang-format on
    return dfeoptions;
}

po::options_description
generate_testing_options(options &o)
{
    po::options_description testing("Testing options");
    // clang-format off
    testing.add_options()("leaf_test",po::bool_switch(&o.leaf_test),"Perform very expensive checking on sample list ranges vs. leaf counts")
        ("matrix_test",po::bool_switch(&o.matrix_test),"Perform run-time test on generating fwdpp::data_matrix objects and validating the row sums")
		("serialization_test",po::value<std::string>(&o.filename),"Test round-trip to/from a file");
    // clang-format on
    return testing;
}

void
validate_primary_options(const options &o)
{
    if (o.theta < 0. || o.rho < 0.)
        {
            throw std::invalid_argument("rho and theta must be >= 0.0");
        }
    if (o.N < 1)
        {
            throw std::invalid_argument("N must be > 0");
        }
    if (o.gcint < 1)
        {
            throw std::invalid_argument(
                "Simplification (gc) interval must be > 0");
        }
    if (o.mu < 0)
        {
            throw std::invalid_argument(
                "Mutation rate to selected variants must be >= 0");
        }
    else if (o.mu > 0)
        {
            if (o.mean == 0.0 && std::isnan(o.scoeff))
                {
                    throw std::invalid_argument(
                        "mean selection coefficient cannot be zero");
                }
        }
    if (o.ancient_sampling_interval > 0 && o.ancient_sample_size < 1)
        {
            throw std::invalid_argument(
                "ansam must be > 0 when tracking ancient samples");
        }
}

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
                                              true, false);
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
                                                        "loopback "
                                                        "error");
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

void
execute_expensive_leaf_test(const options &o,
                            const fwdpp::ts::table_collection &tables,
                            const std::vector<fwdpp::ts::TS_NODE_INT> &samples)
{
    if (o.leaf_test)
        {
            std::cerr << "Starting sample list validation.  This may take a "
                         "while!\n";
            expensive_leaf_test(tables, samples);
            std::cout << "Passed with respect to last generation.\n";
            expensive_leaf_test(tables, tables.preserved_nodes);
            std::cout << "Passed with respect to preserved samples.\n";
        }
}

template <typename poptype>
void
execute_matrix_test_detail(const options &o, const poptype &pop,
                           const fwdpp::ts::table_collection &tables,
                           const std::vector<fwdpp::ts::TS_NODE_INT> &samples)
{
    if (o.matrix_test)
        {
            std::cerr << "Matrix test with respect to last generation...";
            matrix_runtime_test(tables, samples, pop.mutations, pop.mcounts);
            std::cerr << "passed\n";
            if (!tables.preserved_nodes.empty())
                {
                    std::cout << "Matrix test with respect to preserved "
                                 "samples...";
                    matrix_runtime_test(tables, tables.preserved_nodes,
                                        pop.mutations,
                                        pop.mcounts_from_preserved_nodes);
                    std::cerr << "passed\n";
                    auto sc = samples;
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
}

void
execute_matrix_test(const options &o, const single_locus_poptype &pop,
                    const fwdpp::ts::table_collection &tables,
                    const std::vector<fwdpp::ts::TS_NODE_INT> &samples)
{
    execute_matrix_test_detail(o, pop, tables, samples);
}


void
execute_serialization_test(const options &o,
                           const fwdpp::ts::table_collection &tables)
{
    if (!o.filename.empty())
        {
            test_serialization(tables, o.filename);
        }
}

void
write_sfs(const options &o, const fwdpp::GSLrng_mt &rng,
          const fwdpp::ts::table_collection &tables,
          const std::vector<fwdpp::ts::TS_NODE_INT> &samples,
          const std::vector<fwdpp::popgenmut> &mutations)
{
    if (!o.sfsfilename.empty())
        {
            if (!(o.nsam > 2))
                {
                    throw std::invalid_argument(
                        "sample size for site frequency spectrum must be "
                        "> 2");
                }
            // Simplify w.r.to 100 samples
            std::vector<fwdpp::ts::TS_NODE_INT> small_sample(o.nsam);
            auto s(samples);
            gsl_ran_choose(rng.get(), small_sample.data(), small_sample.size(),
                           s.data(), s.size(), sizeof(fwdpp::ts::TS_NODE_INT));
            std::iota(small_sample.begin(), small_sample.end(), 0);
            auto dm = fwdpp::ts::generate_data_matrix(tables, small_sample,
                                                      mutations, true, false, true);
            auto rs = fwdpp::row_sums(dm);
            std::vector<int> sfs(small_sample.size() - 1);
            for (auto i : rs.first)
                {
                    sfs[i - 1]++;
                }
            std::ofstream sfs_stream(o.sfsfilename.c_str());
            for (std::size_t i = 0; i < sfs.size(); ++i)
                {
                    sfs_stream << (i + 1) << ' ' << sfs[i] << '\n';
                }
        }
}
