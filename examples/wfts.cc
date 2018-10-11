/*! \include wfts.cet
 */

#include <cstdio>
#include <iostream>
#include <cassert>
#include <fwdpp/ts/table_collection.hpp>
#include <fwdpp/ts/table_simplifier.hpp>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/sugar/popgenmut.hpp>
#include <fwdpp/sugar/slocuspop.hpp>
#include <fwdpp/fitness_models.hpp>
#include <fwdpp/extensions/callbacks.hpp>
#include <fwdpp/poisson_xover.hpp>
#include <fwdpp/recbinder.hpp>
#include <boost/program_options.hpp>

#include "evolve_generation_ts.hpp"

namespace po = boost::program_options;
using poptype = fwdpp::slocuspop<fwdpp::popgenmut>;
using GSLrng = fwdpp::GSLrng_t<fwdpp::GSL_RNG_MT19937>;

inline fwdpp::fwdpp_internal::gsl_ran_discrete_t_ptr
mean_fitness_zero_out_gametes(poptype &pop)
{
    auto N_curr = pop.diploids.size();
    std::vector<double> fitnesses(N_curr);
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
void
confirm_mutation_counts(poptype &pop,
                        const fwdpp::ts::table_collection &tables,
                        const unsigned generation)
{
#ifndef NDEBUG
    std::vector<std::size_t> keys;
    for (auto &mr : tables.mutation_table)
        {
            keys.push_back(mr.key);
        }
    std::sort(keys.begin(), keys.end());
    auto u = std::unique(keys.begin(), keys.end());
    ;
    if (u != keys.end())
        {
            std::cout << "redundant keys " << generation << '\n';
        }

    for (auto &mr : tables.mutation_table)
        {
            assert(mr.node != fwdpp::ts::TS_NULL_NODE);
            assert(tables.node_table[mr.node].generation
                   >= pop.mutations[mr.key].g);
        }
    decltype(pop.mcounts) mc;
    fwdpp::fwdpp_internal::process_gametes(pop.gametes, pop.mutations, mc);
    for (std::size_t i = 0; i < pop.mcounts.size(); ++i)
        {
            if (!pop.mutations[i].neutral)
                {
                    assert(pop.mcounts[i] == mc[i]);
                }
        }
#endif
}

template <typename poptype>
void
simplify_tables(poptype &pop, fwdpp::ts::table_collection &tables,
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
    tables.count_mutations(pop.mutations, samples, pop.mcounts);
    // TODO: the following steps all need to be updated
    // to deal with mutation counts due to preserved nodes.
    tables.mutation_table.erase(
        std::remove_if(
            tables.mutation_table.begin(), tables.mutation_table.end(),
            [&pop](const fwdpp::ts::mutation_record &mr) {
                return pop.mcounts[mr.key] == 2 * pop.diploids.size();
            }),
        tables.mutation_table.end());
    fwdpp::fwdpp_internal::gamete_cleaner(pop.gametes, pop.mutations,
                                          pop.mcounts, 2 * pop.diploids.size(),
                                          std::true_type());
    fwdpp::update_mutations(pop.mutations, pop.fixations, pop.fixation_times,
                            pop.mut_lookup, pop.mcounts, generation,
                            2 * pop.diploids.size());
    confirm_mutation_counts(pop, tables, generation);
}

int
main(int argc, char **argv)
{
    fwdpp::uint_t N, gcint = 100;
    double theta, rho, mean, shape, mu;
    unsigned seed = 42;
    po::options_description options("Usage");
    // clang-format off
    options.add_options()("help", "Display help")
        ("N", po::value<unsigned>(&N), "Diploid population size")
        ("gc", po::value<unsigned>(&gcint),
        "Simplification interval. Default is 100 generations.")
        ("theta", po::value<double>(&theta), "4Nu")
        ("rho", po::value<double>(&rho), "4Nr")
        ("mu", po::value<double>(&mu), "mutation rate to selected variants")
        ("mean", po::value<double>(&mean), "Mean 2Ns of Gamma distribution of selection coefficients")
        ("shape", po::value<double>(&shape), "Shape of Gamma distribution of selection coefficients")
        ("seed", po::value<unsigned>(&seed), "Random number seed. Default is 42");
    // clang-format on
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, options), vm);
    po::notify(vm);

    if (vm.count("help"))
        {
            std::cout << options << '\n';
            std::exit(1);
        }

    GSLrng rng(seed);

    poptype pop(2 * N);
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
    for (; generation <= 20 * N; ++generation)
        {
            auto lookup = mean_fitness_zero_out_gametes(pop);
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
                    simplify_tables(pop, tables, simplifier,
                                    tables.num_nodes() - 2 * N, 2 * N,
                                    generation);
                    mutation_recycling_bin
                        = fwdpp::fwdpp_internal::make_mut_queue(pop.mcounts);
                    simplified = true;
                    next_index = tables.num_nodes();
                    first_parental_index = 0;
                    confirm_mutation_counts(pop, tables, generation);
                }
            else
                {
                    simplified = false;
                    first_parental_index = next_index;
                    next_index += 2 * N;
                    mutation_recycling_bin = {};
                }
        }
    if (!simplified)
        {
            simplify_tables(pop, tables, simplifier,
                            tables.num_nodes() - 2 * N, 2 * N, generation);
            confirm_mutation_counts(pop, tables, generation);
        }
}
