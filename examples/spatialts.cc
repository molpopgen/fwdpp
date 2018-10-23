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
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/sugar/popgenmut.hpp>
#include <fwdpp/sugar/slocuspop.hpp>
#include <fwdpp/util.hpp>
#include <fwdpp/fitness_models.hpp>
#include <fwdpp/extensions/callbacks.hpp>
#include <fwdpp/poisson_xover.hpp>
#include <fwdpp/recbinder.hpp>
#include <boost/program_options.hpp>

// TODO check for these in configure script
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/core/cs.hpp>
#include <boost/geometry/index/rtree.hpp>
// Boost.Range
#include <boost/range.hpp>
// adaptors
#include <boost/range/adaptor/indexed.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include "evolve_generation_ts.hpp"
#include "confirm_mutation_counts.hpp"

namespace po = boost::program_options;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
using poptype = fwdpp::slocuspop<fwdpp::popgenmut>;
using GSLrng = fwdpp::GSLrng_t<fwdpp::GSL_RNG_MT19937>;
using point = bg::model::point<double, 2, boost::geometry::cs::cartesian>;
using point_to_diploid = std::pair<point, std::size_t>;
using rtree_type = bgi::rtree<point_to_diploid, bgi::quadratic<16>>;

inline fwdpp::fwdpp_internal::gsl_ran_discrete_t_ptr
calculate_fitnesses(poptype &pop, std::vector<double> &fitnesses)
{
    auto N_curr = pop.diploids.size();
    fitnesses.resize(N_curr);
    for (size_t i = 0; i < N_curr; ++i)
        {
            fitnesses[i] = fwdpp::multiplicative_diploid(2.0)(
                pop.diploids[i], pop.gametes, pop.mutations);
        }
    auto lookup = fwdpp::fwdpp_internal::gsl_ran_discrete_t_ptr(
        gsl_ran_discrete_preproc(N_curr, &fitnesses[0]));
    return lookup;
}

template <typename poptype>
std::vector<fwdpp::ts::TS_NODE_INT>
simplify_tables(poptype &pop,
                std::vector<fwdpp::uint_t> &mcounts_from_preserved_nodes,
                fwdpp::ts::table_collection &tables,
                fwdpp::ts::table_simplifier &simplifier,
                const fwdpp::ts::TS_NODE_INT first_sample_node,
                const std::size_t num_samples, const unsigned generation)
{
    tables.sort_tables(pop.mutations);
    std::vector<std::int32_t> samples(num_samples);
    std::iota(samples.begin(), samples.end(), first_sample_node);
    auto idmap = simplifier.simplify(tables, samples, pop.mutations);
    tables.build_indexes();
    for (auto &s : samples)
        {
            s = idmap[s];
        }
    for (auto &s : tables.preserved_nodes)
        {
            assert(idmap[s] != 1);
        }
    fwdpp::ts::count_mutations(tables, pop.mutations, samples, pop.mcounts,
                               mcounts_from_preserved_nodes);
    tables.mutation_table.erase(
        std::remove_if(
            tables.mutation_table.begin(), tables.mutation_table.end(),
            [&pop, &mcounts_from_preserved_nodes](
                const fwdpp::ts::mutation_record &mr) {
                return pop.mcounts[mr.key] == 2 * pop.diploids.size()
                       && mcounts_from_preserved_nodes[mr.key] == 0;
            }),
        tables.mutation_table.end());
    fwdpp::ts::remove_fixations_from_gametes(
        pop.gametes, pop.mutations, pop.mcounts, mcounts_from_preserved_nodes,
        2 * pop.diploids.size());

    fwdpp::ts::flag_mutations_for_recycling(
        pop.mutations, pop.mcounts, mcounts_from_preserved_nodes,
        pop.mut_lookup, 2 * pop.diploids.size());
    confirm_mutation_counts(pop, tables);
    return idmap;
}

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
// From
// https://www.boost.org/doc/libs/1_64_0/libs/geometry/doc/html/geometry/spatial_indexes/rtree_examples/range_adaptors.htmlw
// Define a function object converting a value_type of indexed Range into std::pair<>.
// This is a generic implementation but of course it'd be possible to use some
// specific types. One could also take Value as template parameter and access
// first_type and second_type members, etc.
template <typename First, typename Second> struct pair_maker
{
    typedef std::pair<First, Second> result_type;
    template <typename T>
    inline result_type
    operator()(T const &v) const
    {
        return result_type(v.value(), v.index());
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
    bool leaf_test = false;
    bool matrix_test = false;
    std::string filename, sfsfilename;
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
    // distributed in a circle from (0.5,0.5).  That circle has a
    // radius of 0.5
    std::vector<point> parental_points, offspring_points(N);
    for (fwdpp::uint_t i = 0; i < N; ++i)
        {
            double x, y;
            // This function returns x and y
            // w.r.to a normalized 2d surface
            // such that x^2 + y^2 = 1.
            // Thus, dividing x and y
            // by sqrt(2) gives us x^2+y^2 = 0.5
            gsl_ran_dir_2d(rng.get(), &x, &y);
            x /= std::sqrt(2);
            y /= std::sqrt(2);
            point p{ 0.5 + x, 0.5 + y };
            parental_points.emplace_back(p);
        }
    rtree_type rtree(
        parental_points | boost::adaptors::indexed()
        | boost::adaptors::transformed(pair_maker<point, std::size_t>()));
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
    std::vector<fwdpp::uint_t> mcounts_from_preserved_nodes;
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
              double x = (bg::get<0>(parental_points[p1])
                          + bg::get<0>(parental_points[p2]))
                         / 2.0;
              double y = (bg::get<1>(parental_points[p1])
                          + bg::get<1>(parental_points[p2]))
                         / 2.0;
              // Don't fall of the map...
              if (x < 0.0)
                  {
                      x = 0.0;
                  }
              if (x > 1.0)
                  {
                      x = 1.0;
                  }
              if (y < 0.0)
                  {
                      y = 0.0;
                  }
              if (y > 1.0)
                  {
                      y = 1.0;
                  }
              offspring_points[o] = point(x, y);
          };
    for (; generation <= 10 * N; ++generation)
        {
            //Clear out offspring coordinates
            offspring_points.resize(N);

            auto lookup = calculate_fitnesses(pop, fitnesses);
            auto pick1 = [&lookup, &rng]() {
                return gsl_ran_discrete(rng.get(), lookup.get());
            };
            // NOTE: pick2 is the performance killer!
            auto pick2 = [&rng, &rtree, &parental_points, &fitnesses,
                          mating_radius](const std::size_t p1) {
                std::vector<point_to_diploid> possible_mates;
                //find all individuals in population whose Euclidiean distance
                //from parent1 is <= radius.  The "point" info fill up
                //the possible_mates vector.
                auto pp = parental_points.at(p1);
                rtree.query(bgi::satisfies([pp, mating_radius](
                                               const point_to_diploid &v) {
                                double p1x = bg::get<0>(pp);
                                double p1y = bg::get<1>(pp);
                                double p2x = bg::get<0>(v.first);
                                double p2y = bg::get<1>(v.first);
                                double euclid
                                    = std::sqrt(std::pow(p1x - p2x, 2.0)
                                                + std::pow(p1y - p2y, 2.0));
                                return euclid <= mating_radius;
                            }),
                            std::back_inserter(possible_mates));
                if (possible_mates.size() == 1)
                    return p1; //only possible mate was itself, so we self-fertilize
                // Use inverse CDF to choose parent2 proportional
                // to fitness w/in the mating circle
                // TODO: this seems more cumbersome than necessary
                std::vector<double> fitness_in_radius;
                for (auto pm : possible_mates)
                    {
                        fitness_in_radius.push_back(fitnesses[pm.second]);
                    }
                double sumw = std::accumulate(fitness_in_radius.begin(),
                                              fitness_in_radius.end(), 0.);
                std::vector<double> cumw(fitness_in_radius.size());
                std::partial_sum(fitness_in_radius.begin(),
                                 fitness_in_radius.end(), cumw.begin());
                std::transform(cumw.begin(), cumw.end(), cumw.begin(),
                               [sumw](double w) { return w / sumw; });
                auto x = gsl_rng_uniform(rng.get());
                std::size_t i = 0;
                for (; i < cumw.size(); ++i)
                    {
                        if (cumw[i] > x)
                            break;
                    }
                return possible_mates[i].second;
            };
            evolve_generation(rng, pop, N, mu, pick1, pick2, update_offspring,
                              mmodel, mutation_recycling_bin, recmap,
                              generation, tables, first_parental_index,
                              next_index);
            //rebuild our rtree
            //rtree.clear();
            //for (std::size_t i = 0; i < pop.diploids.size(); ++i)
            //    {
            //        rtree.insert(point_to_diploid{ offspring_points[i], i });
            //    }

            rtree = rtree_type(offspring_points | boost::adaptors::indexed()
                               | boost::adaptors::transformed(
                                     pair_maker<point, std::size_t>()));
            if (generation % gcint == 0.0)
                {
                    std::cout << generation << '\n';
                    auto idmap = simplify_tables(
                        pop, mcounts_from_preserved_nodes, tables, simplifier,
                        tables.num_nodes() - 2 * N, 2 * N, generation);
                    mutation_recycling_bin = fwdpp::ts::make_mut_queue(
                        pop.mcounts, mcounts_from_preserved_nodes);
                    simplified = true;
                    next_index = tables.num_nodes();
                    first_parental_index = 0;
                    confirm_mutation_counts(pop, tables);

                    // When tracking ancient samples, the node ids of those samples change.
                    // Thus, we need to remap our metadata upon simplification
                    for (auto &md : ancient_sample_metadata)
                        {
                            md.n1 = idmap[md.n1];
                            md.n2 = idmap[md.n2];
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
                                i, generation, bg::get<0>(offspring_points[i]),
                                bg::get<1>(offspring_points[i]), x.first,
                                x.second);
                        }
                }
            parental_points.swap(offspring_points);
        }
    if (!simplified)
        {
            auto idmap = simplify_tables(
                pop, mcounts_from_preserved_nodes, tables, simplifier,
                tables.num_nodes() - 2 * N, 2 * N, generation);
            confirm_mutation_counts(pop, tables);
            // When tracking ancient samples, the node ids of those samples change.
            // Thus, we need to remap our metadata upon simplification
            for (auto &md : ancient_sample_metadata)
                {
                    md.n1 = idmap[md.n1];
                    md.n2 = idmap[md.n2];
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
                               mcounts_from_preserved_nodes);
    for (std::size_t i = 0; i < pop.mutations.size(); ++i)
        {
            if (pop.mutations[i].neutral)
                {
                    if (!pop.mcounts[i] && !mcounts_from_preserved_nodes[i])
                        {
                            throw std::runtime_error(
                                "invalid final mutation count");
                        }
                }
        }
}
