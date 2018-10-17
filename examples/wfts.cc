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
#include <fwdpp/ts/table_collection.hpp>
#include <fwdpp/ts/table_simplifier.hpp>
#include <fwdpp/ts/count_mutations.hpp>
#include <fwdpp/ts/recycling.hpp>
#include <fwdpp/ts/marginal_tree_iterator.hpp>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/sugar/popgenmut.hpp>
#include <fwdpp/sugar/slocuspop.hpp>
#include <fwdpp/util.hpp>
#include <fwdpp/fitness_models.hpp>
#include <fwdpp/extensions/callbacks.hpp>
#include <fwdpp/poisson_xover.hpp>
#include <fwdpp/recbinder.hpp>
#include <boost/program_options.hpp>

#include "evolve_generation_ts.hpp"
#include "confirm_mutation_counts.hpp"

namespace po = boost::program_options;
using poptype = fwdpp::slocuspop<fwdpp::popgenmut>;
using GSLrng = fwdpp::GSLrng_t<fwdpp::GSL_RNG_MT19937>;

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

template <typename gcont_t, typename mcont_t,
          typename mutation_count_container>
void
gamete_cleaner(gcont_t &gametes, const mcont_t &mutations,
               const mutation_count_container &mcounts,
               const mutation_count_container &mcounts_from_preserved_nodes,
               const fwdpp::uint_t twoN)
{
    bool fixations_exist = false;
    for (std::size_t i = 0; !fixations_exist && i < mcounts.size(); ++i)
        {
            if (mcounts[i] == twoN && mcounts_from_preserved_nodes[i] == 0)
                {
                    fixations_exist = true;
                }
        }
    if (fixations_exist)
        {
            auto removal_criteria = [&mcounts, &mcounts_from_preserved_nodes,
                                     twoN](const fwdpp::uint_t key) {
                return mcounts[key] == twoN
                       && mcounts_from_preserved_nodes[key] == 0;
            };
            for (auto &g : gametes)
                {
                    if (g.n)
                        {
                            auto itr = std::remove_if(g.mutations.begin(),
                                                      g.mutations.end(),
                                                      removal_criteria);
                            g.mutations.erase(itr, g.mutations.end());
                            itr = std::remove_if(g.smutations.begin(),
                                                 g.smutations.end(),
                                                 removal_criteria);
                            g.smutations.erase(itr, g.smutations.end());
                        }
                }
        }
}

template <typename mcont_t, typename lookup_table,
          typename mutation_count_container>
void
update_mutations(const mcont_t &mutations, mutation_count_container &mcounts,
                 mutation_count_container &mcounts_from_preserved_nodes,
                 lookup_table &lookup, const fwdpp::uint_t twoN)
{
    for (std::size_t i = 0; i < mcounts.size(); ++i)
        {
            if (mcounts_from_preserved_nodes[i] == 0)
                {
                    if (mcounts[i] > twoN)
                        {
                            throw std::runtime_error(
                                "mutation count out of range");
                        }
                    if (mcounts[i] == twoN)
                        {
                            auto itr = lookup.equal_range(mutations[i].pos);
                            while (itr.first != itr.second)
                                {
                                    if (itr.first->second == i)
                                        {
                                            lookup.erase(itr.first);
                                            mcounts[i] = 0;
                                            break;
                                        }
                                    ++itr.first;
                                }
                        }
                    else if (mcounts[i] == 0)
                        {
                            auto itr = lookup.equal_range(mutations[i].pos);
                            if (itr.first != lookup.end())
                                {
                                    while (itr.first != itr.second)
                                        {
                                            if (itr.first->second == i)
                                                {
                                                    lookup.erase(itr.first);
                                                    break;
                                                }
                                            ++itr.first;
                                        }
                                }
                        }
                }
        }
}

// TODO: consider flattening the return value to a vector
std::map<fwdpp::ts::TS_NODE_INT, std::vector<std::pair<double, double>>>
mark_multiple_roots(const fwdpp::ts::table_collection &tables,
                    const std::vector<fwdpp::ts::TS_NODE_INT> &samples)
{
    std::map<fwdpp::ts::TS_NODE_INT, std::vector<std::pair<double, double>>>
        rv;
    fwdpp::ts::marginal_tree_iterator mti(tables, samples);
    while (mti(std::true_type(), std::false_type()))
        {
            bool single_root = false;
            for (auto &s : samples)
                {
                    auto p = s;
                    auto lp = p;
                    while (p != -1)
                        {
                            lp = p;
                            p = mti.marginal.parents[p];
                        }
                    if (mti.marginal.leaf_counts[lp] == samples.size())
                        {
                            single_root = true;
                        }
                    else
                        {
                            auto itr = rv.find(lp);
                            auto w = std::make_pair(mti.marginal.left,
                                                    mti.marginal.right);
                            if (itr == rv.end())
                                {
                                    rv[lp].emplace_back(std::move(w));
                                }
                            else
                                {
                                    assert(!itr->second.empty());
                                    if (std::find(itr->second.begin(),
                                                  itr->second.end(), w)
                                        == itr->second.end())
                                        {
                                            itr->second.emplace_back(
                                                std::move(w));
                                        }
                                }
                        }
                    if (single_root)
                        {
                            break;
                        }
                }
        }
    return rv;
}

template <typename rng, typename mfunction>
unsigned
mutate_tables(const rng &r, const mfunction &make_mutation,
              fwdpp::ts::table_collection &tables,
              const std::vector<fwdpp::ts::TS_NODE_INT> &samples,
              const double mu)
{
    unsigned nmuts = 0;
    auto mr = mark_multiple_roots(tables, samples);
    for (auto &i : mr)
        {
            auto dt = tables.node_table[i.first].generation;
            for (auto j : i.second)
                {
                    double mean = dt * (j.second - j.first) * mu;
                    auto nm = gsl_ran_poisson(r.get(), mean);
                    nmuts += nm;
                    for (unsigned m = 0; m < nm; ++m)
                        {
                            unsigned g = static_cast<unsigned>(
                                gsl_ran_flat(r.get(), 1, dt + 1));
                            auto k = make_mutation(j.first, j.second, g);
                            tables.mutation_table.emplace_back(
                                fwdpp::ts::mutation_record{ i.first, k });
                        }
                }
        }
    for (auto &e : tables.edge_table)
        {
            auto ct = tables.node_table[e.child].generation;
            auto pt = tables.node_table[e.parent].generation;
            auto dt = ct - pt;
            double mean = dt * (e.right - e.left) * mu;
            auto nm = gsl_ran_poisson(r.get(), mean);
            for (unsigned m = 0; m < nm; ++m)
                {
                    unsigned g = static_cast<unsigned>(
                        gsl_ran_flat(r.get(), pt + 1, ct + 1));
                    auto k = make_mutation(e.left, e.right, g);
                    tables.mutation_table.emplace_back(
                        fwdpp::ts::mutation_record{ e.child, k });
                }
            nmuts += nm;
        }
    return nmuts;
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
    gamete_cleaner(pop.gametes, pop.mutations, pop.mcounts,
                   mcounts_from_preserved_nodes, 2 * pop.diploids.size());

    update_mutations(pop.mutations, pop.mcounts, mcounts_from_preserved_nodes,
                     pop.mut_lookup, 2 * pop.diploids.size());
    confirm_mutation_counts(pop, tables);
    return idmap;
}

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

int
main(int argc, char **argv)
{
    fwdpp::uint_t N, gcint = 100;
    double theta, rho, mean = 0.0, shape = 1, mu;
    unsigned seed = 42;
    int ancient_sampling_interval = -1;
    int ancient_sample_size = -1;
    po::options_description options("Usage");
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
         "Sample size (no. diploids) of ancient samples to take at each ancient sampling interval.  Default is -1, and must be reset if sampling_interval is used");
    // clang-format on
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, options), vm);
    po::notify(vm);

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
    if (vm.count("help"))
        {
            std::cout << options << '\n';
            std::exit(1);
        }

    GSLrng rng(seed);

    poptype pop(N);
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
    for (; generation <= 10 * N; ++generation)
        {
            auto lookup = calculate_fitnesses(pop, fitnesses);
            auto pick1 = [&lookup, &rng]() {
                return gsl_ran_discrete(rng.get(), lookup.get());
            };
            auto pick2 = [&lookup, &rng](const std::size_t /*p1*/) {
                return gsl_ran_discrete(rng.get(), lookup.get());
            };
            evolve_generation(rng, pop, N, mu, pick1, pick2, mmodel,
                              mutation_recycling_bin, recmap, generation,
                              tables, first_parental_index, next_index);
            if (generation % gcint == 0.0)
                {
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
                            assert(tables.node_table[x.first].generation
                                   == generation);
                            assert(tables.node_table[x.second].generation
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
                                i, generation, fitnesses[i], x.first,
                                x.second);
                        }
                }
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
            if (tables.node_table[mr.n1].generation != mr.time
                || tables.node_table[mr.n2].generation != mr.time)
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
    auto neutral_muts = mutate_tables(rng, neutral_variant_maker, tables, s,
                                      theta / static_cast<double>(4 * N));
    std::sort(tables.mutation_table.begin(), tables.mutation_table.end(),
              [&pop](const fwdpp::ts::mutation_record &a,
                     const fwdpp::ts::mutation_record &b) {
                  return pop.mutations[a.key].pos < pop.mutations[b.key].pos;
              });
    fwdpp::ts::count_mutations(tables, pop.mutations, s, pop.mcounts);
    for (std::size_t i = 0; i < pop.mutations.size(); ++i)
        {
            if (pop.mutations[i].neutral)
                {
                    if (!pop.mcounts[i])
                        {
                            throw std::runtime_error(
                                "invalid final mutation count");
                        }
                }
        }
    std::cout << neutral_muts << '\n';
}
