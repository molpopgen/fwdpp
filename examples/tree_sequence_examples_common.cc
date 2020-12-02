#include <cmath>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <gsl/gsl_randist.h>
#include <fwdpp/ts/count_mutations.hpp>
#include <fwdpp/ts/generate_data_matrix.hpp>
#include <fwdpp/ts/serialization.hpp>
#include <fwdpp/ts/mutate_tables.hpp>
#include <fwdpp/ts/mutation_tools.hpp>
#include <fwdpp/ts/marginal_tree_functions/samples.hpp>
#include <fwdpp/ts/visit_sites.hpp>
#include <fwdpp/extensions/callbacks.hpp>
#include "tree_sequence_examples_common.hpp"

namespace po = boost::program_options;

template <typename poptype>
int
apply_neutral_mutations_details(const options &o, const GSLrng &rng,
                                fwdpp::ts::std_table_collection &tables, poptype &pop,
                                fwdpp::flagged_mutation_queue &mutation_recycling_bin)
{
    std::vector<fwdpp::ts::table_index_t> s(2 * o.N);

    std::iota(s.begin(), s.end(), 0);
    const auto neutral_variant_maker = [&rng, &pop, &mutation_recycling_bin](
                                           const double left, const double right,
                                           const fwdpp::uint_t generation) {
        auto key = fwdpp::infsites_mutation(
            mutation_recycling_bin, pop.mutations, rng.get(), pop.mut_lookup, generation,
            0.0, [left, right, &rng] { return gsl_ran_flat(rng.get(), left, right); },
            []() { return 0.0; }, []() { return 0.0; });
        return fwdpp::ts::new_variant_record(pop.mutations[key].pos, 0, key,
                                             pop.mutations[key].neutral, 1);
    };
    auto neutral_muts = fwdpp::ts::mutate_tables(rng, neutral_variant_maker, tables, s,
                                                 o.theta / static_cast<double>(4 * o.N));
    return neutral_muts;
}

int
apply_neutral_mutations(const options &o, const GSLrng &rng,
                        fwdpp::ts::std_table_collection &tables,
                        ts_examples_poptype &pop,
                        fwdpp::flagged_mutation_queue &mutation_recycling_bin)
{
    return apply_neutral_mutations_details(o, rng, tables, pop, mutation_recycling_bin);
}

options::options()
    : N{}, gcint(100), theta(), rho(), mean(0.), shape(1.), mu(),
      scoeff(std::numeric_limits<double>::quiet_NaN()), dominance(1.), scaling(2.),
      seed(42), ancient_sampling_interval(-1), ancient_sample_size(-1), nsam(0),
      leaf_test(false), matrix_test(false), visit_sites_test(false),
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
		("serialization_test",po::value<std::string>(&o.filename),"Test round-trip to/from a file")
		("visit_sites_tests",po::bool_switch(&o.visit_sites_test),"Test correctness of ts::visit_sites");
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
            throw std::invalid_argument("Simplification (gc) interval must be > 0");
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
make_dfe(const fwdpp::uint_t N, const GSLrng &r, const double mean, const double shape,
         const double scoeff)
{
    if (std::isfinite(scoeff))
        {
            return [scoeff]() { return scoeff; };
        }
    fwdpp::extensions::gamma dfe(mean, shape);
    return [&r, dfe, N]() { return dfe(r.get()) / static_cast<double>(2 * N); };
}

void
matrix_runtime_test(const fwdpp::ts::std_table_collection &tables,
                    const std::vector<fwdpp::ts::table_index_t> &samples,
                    const std::vector<fwdpp::uint_t> &mcounts)
{
    auto dm = fwdpp::ts::generate_data_matrix(tables, samples, true, true, false);
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
expensive_leaf_test(const fwdpp::ts::std_table_collection &tables,
                    const std::vector<fwdpp::ts::table_index_t> &sample_list)
{
    fwdpp::ts::tree_visitor<fwdpp::ts::std_table_collection> mti(
        tables, sample_list, fwdpp::ts::update_samples_list(true));
    while (mti())
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
                                                    throw std::runtime_error("loopback "
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
test_serialization(const fwdpp::ts::std_table_collection &tables,
                   const std::string &filename)
{
    std::ofstream o(filename.c_str());
    fwdpp::ts::io::serialize_tables(o, tables);
    o.close();
    std::ifstream i(filename.c_str());
    auto tables2
        = fwdpp::ts::io::deserialize_tables<fwdpp::ts::std_table_collection>()(i);

    if (tables.genome_length() != tables2.genome_length())
        {
            throw std::runtime_error("genome_length does not match");
        }
    if (tables.edge_offset != tables2.edge_offset)
        {
            throw std::runtime_error("edge_offset does not match");
        }

    if (tables.edges != tables2.edges)
        {
            throw std::runtime_error("edge tables do not match");
        }

    if (tables.nodes != tables2.nodes)
        {
            throw std::runtime_error("node tables do not match");
        }

    if (tables.mutations != tables2.mutations)
        {
            throw std::runtime_error("mutation tables do not match");
        }
    if (tables.sites != tables2.sites)
        {
            throw std::runtime_error("site tables do no match");
        }
    if (tables != tables2)
        {
            throw std::runtime_error("tables failed equality check");
        }
    std::cout << "serialization test passed\n";
}

void
execute_expensive_leaf_test(const options &o,
                            const fwdpp::ts::std_table_collection &tables,
                            const std::vector<fwdpp::ts::table_index_t> &samples,
                            const std::vector<fwdpp::ts::table_index_t> &preserved_nodes)
{
    if (o.leaf_test)
        {
            std::cerr << "Starting sample list validation.  This may take a "
                         "while!\n";
            expensive_leaf_test(tables, samples);
            std::cout << "Passed with respect to last generation.\n";
            expensive_leaf_test(tables, preserved_nodes);
            std::cout << "Passed with respect to preserved samples.\n";
        }
}

template <typename poptype>
void
execute_matrix_test_detail(const options &o, const poptype &pop,
                           const fwdpp::ts::std_table_collection &tables,
                           const std::vector<fwdpp::ts::table_index_t> &samples,
                           const std::vector<fwdpp::ts::table_index_t> &preserved_nodes)
{
    if (o.matrix_test)
        {
            std::cerr << "Matrix test with respect to last generation...";
            matrix_runtime_test(tables, samples, pop.mcounts);
            std::cerr << "passed\n";
            // NOTE: the following check is useful, but disabled by default.
            //decltype(pop.mcounts) mcounts_from_genomes(pop.mutations.size(), 0);
            //for (auto &d : pop.diploids)
            //    {
            //        for (auto k : pop.haploid_genomes[d.first].smutations)
            //            {
            //                mcounts_from_genomes[k]++;
            //            }
            //        for (auto k : pop.haploid_genomes[d.second].smutations)
            //            {
            //                mcounts_from_genomes[k]++;
            //            }
            //    }
            //for (std::size_t i = 0; i < pop.mutations.size(); ++i)
            //    {
            //        if (pop.mutations[i].neutral == false)
            //            {
            //                if (pop.mcounts[i] != mcounts_from_genomes[i])
            //                    {
            //                        throw std::runtime_error(
            //                            "Matrix test: tree sequence mutation count "
            //                            "ddoesn't match data in genomes");
            //                    }
            //            }
            //    }
            if (!preserved_nodes.empty())
                {
                    std::cout << "Matrix test with respect to preserved "
                                 "samples...";
                    matrix_runtime_test(tables, preserved_nodes,
                                        pop.mcounts_from_preserved_nodes);
                    std::cerr << "passed\n";
                    auto sc = samples;
                    sc.insert(sc.end(), preserved_nodes.begin(), preserved_nodes.end());
                    auto mc(pop.mcounts);
                    std::transform(mc.begin(), mc.end(),
                                   pop.mcounts_from_preserved_nodes.begin(), mc.begin(),
                                   std::plus<fwdpp::uint_t>());
                    std::cout << "Matrix test with respect to last generation "
                                 "+ preserved nodes...";
                    matrix_runtime_test(tables, sc, mc);
                    std::cout << "passed.\n";
                    std::cout << "Matrix test with respect to most recent "
                                 "ancient sampling time point...";
                    sc.clear();
                    std::copy_if(
                        preserved_nodes.begin(), preserved_nodes.end(),
                        std::back_inserter(sc),
                        [&tables, &preserved_nodes](const fwdpp::ts::table_index_t n) {
                            return tables.nodes[n].time
                                   == tables.nodes[preserved_nodes.back()].time;
                        });
                    mc.clear();
                    fwdpp::ts::count_mutations(tables, pop.mutations, sc, mc);
                    matrix_runtime_test(tables, sc, mc);
                    std::cout << "passed\n";
                }
        }
}

void
visit_sites_test(const options &o, const ts_examples_poptype &pop,
                 const fwdpp::ts::std_table_collection &tables,
                 const std::vector<fwdpp::ts::table_index_t> &samples,
                 const std::vector<fwdpp::ts::table_index_t> &preserved_nodes)
{
    if (o.visit_sites_test)
        {
            auto mc(pop.mcounts);
            mc.clear();
            auto s(samples);
            s.insert(end(s), begin(preserved_nodes), end(preserved_nodes));
            fwdpp::ts::count_mutations(tables, pop.mutations, samples, mc);
            auto mc2(mc);
            std::fill(begin(mc2), end(mc2), 0);
            auto f =
                [&mc2](const fwdpp::ts::marginal_tree &m, const fwdpp::ts::site & /*s*/,
                       std::vector<fwdpp::ts::mutation_record>::const_iterator b,
                       const std::vector<fwdpp::ts::mutation_record>::const_iterator e) {
                    for (; b < e; ++b)
                        {
                            mc2[b->key] = fwdpp::ts::num_samples(m, b->node);
                        }
                };
            fwdpp::ts::visit_sites(tables, samples, f, 0., tables.genome_length());
            if (mc != mc2)
                {
                    throw std::runtime_error("visit_sites_test failed to "
                                             "correctly traverse all sites");
                }
            std::cout << "visit_sites_test passed\n";
        }
}
void
execute_matrix_test(const options &o, const ts_examples_poptype &pop,
                    const fwdpp::ts::std_table_collection &tables,
                    const std::vector<fwdpp::ts::table_index_t> &samples,
                    const std::vector<fwdpp::ts::table_index_t> &preserved_nodes)
{
    execute_matrix_test_detail(o, pop, tables, samples, preserved_nodes);
}

void
execute_serialization_test(const options &o,
                           const fwdpp::ts::std_table_collection &tables)
{
    if (!o.filename.empty())
        {
            test_serialization(tables, o.filename);
        }
}

void
write_sfs(const options &o, const GSLrng &rng,
          const fwdpp::ts::std_table_collection &tables,
          const std::vector<fwdpp::ts::table_index_t> &samples)
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
            std::vector<fwdpp::ts::table_index_t> small_sample(o.nsam);
            auto s(samples);
            gsl_ran_choose(rng.get(), small_sample.data(), small_sample.size(), s.data(),
                           s.size(), sizeof(fwdpp::ts::table_index_t));
            std::iota(small_sample.begin(), small_sample.end(), 0);
            auto dm = fwdpp::ts::generate_data_matrix(tables, small_sample, true, false,
                                                      true);
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
