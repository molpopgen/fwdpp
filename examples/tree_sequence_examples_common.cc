#include <cmath>
#include <stdexcept>
#include <fstream>
#include <fwdpp/ts/generate_data_matrix.hpp>
#include <fwdpp/ts/serialization.hpp>
#include <fwdpp/extensions/callbacks.hpp>
#include "tree_sequence_examples_common.hpp"

namespace po = boost::program_options;

options::options()
    : N{}, gcint(100), theta(), rho(), mean(0.), shape(1.), mu(),
      scoeff(std::numeric_limits<double>::max()), dominance(1.), scaling(2.),
      seed(42), ancient_sampling_interval(-1), ancient_sample_size(-1),
      nsam(0), leaf_test(false), matrix_test(false), preserve_fixations(false),
      filename(), sfsfilename()
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
