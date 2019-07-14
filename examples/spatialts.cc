/*! \include spatialts.cc
 * Discrete generations on a continuous landscape using tree sequences.
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
#include <fwdpp/ts/mutation_tools.hpp>
#include <fwdpp/GSLrng_t.hpp>
#include <fwdpp/popgenmut.hpp>
#include <fwdpp/diploid_population.hpp>
#include <fwdpp/util.hpp>
#include <fwdpp/fitness_models.hpp>
#include <fwdpp/extensions/callbacks.hpp>
#include <fwdpp/poisson_xover.hpp>
#include <fwdpp/recbinder.hpp>
#include <fwdpp/simparams.hpp>
#include <boost/program_options.hpp>

#include "simplify_tables.hpp"
#include "evolve_generation_ts.hpp"
#include "calculate_fitnesses.hpp"

namespace po = boost::program_options;
using poptype = fwdpp::diploid_population<fwdpp::popgenmut>;
using GSLrng = fwdpp::GSLrng_t<fwdpp::GSL_RNG_MT19937>;

struct location
{
    double x, y;
    std::size_t individual;
    location()
        : x(std::numeric_limits<double>::quiet_NaN()),
          y(std::numeric_limits<double>::quiet_NaN()),
          individual(std::numeric_limits<std::size_t>::max())
    {
    }
    location(double x_, double y_, double i_) : x(x_), y(y_), individual(i_) {}
};

struct diploid_metadata
{
    std::size_t individual;
    double time, x, y;
    fwdpp::ts::TS_NODE_INT n1, n2;
    diploid_metadata(std::size_t i, double t, double x_, double y_,
                     fwdpp::ts::TS_NODE_INT a, fwdpp::ts::TS_NODE_INT b)
        : individual(i), time(t), x(x_), y(y_), n1(a), n2(b)
    {
    }
};

int
main(int argc, char **argv)
{
    fwdpp::uint_t N, gcint = 100;
    double theta, rho, mean = 0.0, shape = 1, mu, mating_radius = 0.1;
    unsigned seed = 42;
    int ancient_sampling_interval = -1;
    int ancient_sample_size = -1, nsam = 0;
    //bool leaf_test = false;
    //bool matrix_test = false;
    std::string sfsfilename;
    po::options_description options("Simulation options"),
        landscape_options("Landscape options");
    //testing("Testing options"),
    // clang-format off
    options.add_options()("help", "Display help")
        ("N", po::value<unsigned>(&N), "Diploid population size")
        ("gc", po::value<unsigned>(&gcint),
        "Simplification interval. Default is 100 generations.")
        ("theta", po::value<double>(&theta), "4Nu")
        ("rho", po::value<double>(&rho), "4Nr")
        ("mu", po::value<double>(&mu), "mutation rate to selected variants")
        ("mean", po::value<double>(&mean), "Mean 2Ns of Gamma distribution of selection coefficients. Default 0.0.")
        ("shape", po::value<double>(&shape), "Shape of Gamma distribution of selection coefficients. Default = 1.")
        ("seed", po::value<unsigned>(&seed), "Random number seed. Default is 42")
        ("sampling_interval", po::value<int>(&ancient_sampling_interval), 
         "How often to preserve ancient samples.  Default is -1, which means do not preserve any.")
        ("ansam", po::value<int>(&ancient_sample_size),
         "Sample size (no. diploids) of ancient samples to take at each ancient sampling interval.  Default is -1, and must be reset if sampling_interval is used")
		("sfs", po::value<std::string>(&sfsfilename),"Write the site frequency spectrum of a sample to a file")
		("nsam", po::value<int>(&nsam), "Sample size for the site frequency spectrum.  Default is 0.  Change when using --sfs");
        landscape_options.add_options()
        ("mating_radius",po::value<double>(&mating_radius),"Mating radius. Default is 0.1");
        //testing.add_options()("leaf_test",po::bool_switch(&leaf_test),"Perform very expensive checking on sample list ranges vs. leaf counts")
        //("matrix_test",po::bool_switch(&matrix_test),"Perform run-time test on generating fwdpp::data_matrix objects and validating the row sums")
		//("serialization_test",po::value<std::string>(&filename),"Test round-trip to/from a file");
    // clang-format on
    options.add(landscape_options);
    //options.add(testing);
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
            if (mean == 0.0)
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

    // Create our landscape, which will be a square on (0,0) to (1,1)
    // Initially, we assign all individuals to a location uniformly
    // distributed along each axis
    std::vector<location> parental_points, offspring_points(N);
    for (fwdpp::uint_t i = 0; i < N; ++i)
        {
            double x = gsl_rng_uniform(rng.get());
            double y = gsl_rng_uniform(rng.get());
            parental_points.emplace_back(x, y, i);
        }
    fwdpp::ts::table_collection tables(2 * pop.diploids.size(), 0, 0, 1.0);
    fwdpp::ts::table_simplifier simplifier(1.0);
    unsigned generation = 1;
    double recrate = rho / static_cast<double>(4 * N);
    const auto recmap
        = fwdpp::recbinder(fwdpp::poisson_xover(recrate, 0., 1.), rng.get());

    const fwdpp::extensions::gamma dfe(mean, shape);
    const auto get_selection_coefficient = [&rng, dfe, N]() {
        return dfe(rng.get()) / static_cast<double>(2 * N);
    };
    const auto generate_mutation_position
        = [&rng]() { return gsl_rng_uniform(rng.get()); };
    const auto generate_h = []() { return 1.0; };
    const auto make_mutation
        = [&pop, &rng, &generation, generate_mutation_position,
           get_selection_coefficient,
           generate_h](fwdpp::flagged_mutation_queue &recbin,
                       poptype::mcont_t &mutations) {
              return fwdpp::infsites_popgenmut(
                  recbin, mutations, rng.get(), pop.mut_lookup, generation,
                  // 1.0 signifies 100% of mutations will be selected
                  1.0, generate_mutation_position, get_selection_coefficient,
                  generate_h);
          };

    const auto mmodel = [&rng, mu,
                         &make_mutation](fwdpp::flagged_mutation_queue &recbin,
                                         poptype::mcont_t &mutations) {
        std::vector<fwdpp::uint_t> rv;
        unsigned nmuts = gsl_ran_poisson(rng.get(), mu);
        for (unsigned i = 0; i < nmuts; ++i)
            {
                rv.push_back(make_mutation(recbin, mutations));
            }
        std::sort(begin(rv), end(rv),
                  [&mutations](const fwdpp::uint_t a, const fwdpp::uint_t b) {
                      return mutations[a].pos < mutations[b].pos;
                  });
        return rv;
    };

    // Evolve pop for 20N generations
    fwdpp::ts::TS_NODE_INT first_parental_index = 0,
                           next_index = 2 * pop.diploids.size();
    bool simplified = false;
    std::vector<std::size_t> individual_labels(N);
    std::iota(individual_labels.begin(), individual_labels.end(), 0);
    std::vector<std::size_t> individuals;
    if (ancient_sample_size > 0)
        {
            individuals.resize(ancient_sample_size);
        }
    std::vector<double> fitnesses;
    std::vector<diploid_metadata> ancient_sample_metadata;
    // Set offspring location = parental midpoint
    // TODO: add dispersal radius
    const auto update_offspring
        = [&offspring_points, &parental_points](std::size_t o, std::size_t p1,
                                                std::size_t p2) {
              auto &opoint = offspring_points[o];
              opoint.individual = o;
              opoint.x = (parental_points[p1].x + parental_points[p2].x) / 2.0;
              opoint.y = (parental_points[p1].y + parental_points[p2].y) / 2.0;
              // Don't fall of the map...
              if (opoint.x < 0.0)
                  {
                      opoint.x = 0.0;
                  }
              if (opoint.x > 1.0)
                  {
                      opoint.x = 1.0;
                  }
              if (opoint.y < 0.0)
                  {
                      opoint.y = 0.0;
                  }
              if (opoint.y > 1.0)
                  {
                      opoint.y = 1.0;
                  }
          };
    std::vector<std::size_t> possible_mates;
    std::vector<double> possible_mate_fitness;
    std::vector<double> cumw;
    auto ff = fwdpp::multiplicative_diploid(fwdpp::fitness(2.0));
    auto genetics = fwdpp::make_genetic_parameters(
        std::move(ff), std::move(mmodel), std::move(recmap));
    for (; generation <= 10 * N; ++generation)
        {
            //Clear out offspring coordinates
            offspring_points.resize(N);

            auto lookup = calculate_fitnesses(pop, fitnesses, ff);
            auto pick1 = [&lookup, &rng]() {
                return gsl_ran_discrete(rng.get(), lookup.get());
            };
            // NOTE: pick2 is the performance killer!
            auto pick2 = [&rng, &parental_points, &possible_mates, &fitnesses,
                          &possible_mate_fitness, &cumw,
                          mating_radius](const std::size_t p1) {
                possible_mates.clear();
                possible_mate_fitness.clear();
                //find all individuals in population whose Euclidiean distance
                //from parent1 is <= radius.  The "point" info fill up
                //the possible_mates vector.
                auto pp = parental_points[p1];
                for (auto &v : parental_points)
                    {
                        double p1x = pp.x;
                        double p1y = pp.y;
                        double p2x = v.x;
                        double p2y = v.y;
                        double euclid = std::sqrt(std::pow(p1x - p2x, 2.0)
                                                  + std::pow(p1y - p2y, 2.0));
                        if (euclid <= mating_radius)
                            {
                                possible_mates.push_back(v.individual);
                                possible_mate_fitness.push_back(
                                    fitnesses[v.individual]);
                            }
                    }
                if (possible_mates.size() == 1)
                    return p1; //only possible mate was itself, so we self-fertilize
                // Use inverse CDF to choose parent2 proportional
                // to fitness w/in the mating circle
                // TODO: this seems more cumbersome than necessary
                cumw.resize(possible_mate_fitness.size());
                std::partial_sum(possible_mate_fitness.begin(),
                                 possible_mate_fitness.end(), cumw.begin());
                auto sumw = cumw.back();
                std::transform(cumw.begin(), cumw.end(), cumw.begin(),
                               [sumw](double w) { return w / sumw; });
                auto x = gsl_rng_uniform(rng.get());
                std::size_t i = 0;
                for (; i < cumw.size(); ++i)
                    {
                        if (cumw[i] > x)
                            break;
                    }
                return possible_mates[i];
            };
            evolve_generation(rng, pop, genetics, N, pick1, pick2,
                              update_offspring, generation, tables,
                              first_parental_index, next_index);

            if (generation % gcint == 0.0)
                {
                    auto rv = simplify_tables(
                        pop, generation, pop.mcounts_from_preserved_nodes,
                        tables, simplifier, tables.num_nodes() - 2 * N, 2 * N);
                    genetics.mutation_recycling_bin
                        = fwdpp::ts::make_mut_queue(
                            pop.mcounts, pop.mcounts_from_preserved_nodes);
                    simplified = true;
                    next_index = tables.num_nodes();
                    first_parental_index = 0;
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
            else
                {
                    simplified = false;
                    first_parental_index = next_index;
                    next_index += 2 * N;
                }
            if (ancient_sampling_interval > 0
                && generation % ancient_sampling_interval == 0.0
                // This last check forbids us recording the
                // final generation as ancient samples.
                && generation < 10 * N)
                {
                    // For recording the metadata, let's normalize the
                    // fitnesses so that we record what matters in the sim,
                    // which is relative fitness.
                    auto wbar = std::accumulate(fitnesses.begin(),
                                                fitnesses.end(), 0.)
                                / static_cast<double>(N);
                    std::transform(fitnesses.begin(), fitnesses.end(),
                                   fitnesses.begin(),
                                   [wbar](double w) { return w / wbar; });
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
                            ancient_sample_metadata.emplace_back(
                                i, generation, offspring_points[i].x,
                                offspring_points[i].y, x.first, x.second);
                        }
                }
            parental_points.swap(offspring_points);
        }
    if (!simplified)
        {
            auto rv = simplify_tables(
                pop, generation, pop.mcounts_from_preserved_nodes, tables,
                simplifier, tables.num_nodes() - 2 * N, 2 * N);
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
    const auto neutral_variant_maker =
        [&rng, &pop, &genetics](const double left, const double right,
                                const fwdpp::uint_t generation) {
            auto key = fwdpp::infsites_popgenmut(
                genetics.mutation_recycling_bin, pop.mutations, rng.get(),
                pop.mut_lookup, generation, 0.0,
                [left, right, &rng] {
                    return gsl_ran_flat(rng.get(), left, right);
                },
                []() { return 0.0; }, []() { return 0.0; });
            return fwdpp::ts::new_variant_record(
                pop.mutations[key].pos, 0, key, pop.mutations[key].neutral, 1);
        };
    auto neutral_muts
        = fwdpp::ts::mutate_tables(rng, neutral_variant_maker, tables, s,
                                   theta / static_cast<double>(4 * N));
    fwdpp::ts::count_mutations(tables, pop.mutations, s, pop.mcounts,
                               pop.mcounts_from_preserved_nodes);
    for (std::size_t i = 0; i < pop.mutations.size(); ++i)
        {
            if (pop.mutations[i].neutral)
                {
                    if (!pop.mcounts[i]
                        && !pop.mcounts_from_preserved_nodes[i])
                        {
                            throw std::runtime_error(
                                "invalid final mutation count");
                        }
                }
        }
}
