/*! \include wfts.cc
 * Wright-Fisher simulation with tree sequences.
 *
 * NOTE: this program is not streamlined for maximal efficiency.
 * Rather, it serves as a stochastic integration test.  Lots of the
 * inner workings get stress-tested repeatedly during execution.
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
#include <fwdpp/ts/recycling.hpp>
#include <fwdpp/ts/remove_fixations_from_gametes.hpp>
#include <fwdpp/util.hpp>
#include <fwdpp/fitness_models.hpp>
#include <fwdpp/extensions/callbacks.hpp>
#include <fwdpp/genetic_map/genetic_map.hpp>
#include <fwdpp/genetic_map/poisson_interval.hpp>
#include <fwdpp/recbinder.hpp>
#include <fwdpp/simparams.hpp>

#include "simplify_tables.hpp"
#include "evolve_generation_ts.hpp"
#include "calculate_fitnesses.hpp"
#include "tree_sequence_examples_common.hpp"

namespace po = boost::program_options;

int
main(int argc, char **argv)
{
    options o;

    auto main_options = generate_main_options(o);
    auto dfe_options = generate_dfe_options(o);
    auto testing_options = generate_testing_options(o);
    main_options.add(dfe_options);
    main_options.add(testing_options);
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, main_options), vm);
    po::notify(vm);

    if (vm.count("help") || argc == 1)
        {
            std::cout << main_options << '\n';
            std::exit(1);
        }

    validate_primary_options(o);

    GSLrng rng(o.seed);

    ts_examples_poptype pop(o.N);
    fwdpp::ts::std_table_collection tables(2 * pop.diploids.size(), 0, 0, 1.0);
    auto simplifier = fwdpp::ts::make_table_simplifier(tables);
    unsigned generation = 1;
    double recrate = o.rho / static_cast<double>(4 * o.N);
    fwdpp::genetic_map gm;
    gm.add_callback(fwdpp::poisson_interval(0., 1., recrate));
    const auto recmap = fwdpp::recbinder(std::cref(gm), rng.get());

    auto get_selection_coefficient = make_dfe(o.N, rng, o.mean, o.shape, o.scoeff);
    const auto generate_mutation_position
        = [&rng]() { return gsl_rng_uniform(rng.get()); };
    const auto generate_h = [&o]() { return o.dominance; };
    const auto make_mutation =
        [&pop, &rng, &generation, generate_mutation_position, get_selection_coefficient,
         generate_h](fwdpp::flagged_mutation_queue &recbin,
                     ts_examples_poptype::mutation_container &mutations) {
            return fwdpp::infsites_mutation(
                recbin, mutations, rng.get(), pop.mut_lookup, generation,
                // 1.0 signifies 100% of mutations will be selected
                1.0, generate_mutation_position, get_selection_coefficient, generate_h);
        };

    const auto mmodel =
        [&rng, &o, &make_mutation](fwdpp::flagged_mutation_queue &recbin,
                                   ts_examples_poptype::mutation_container &mutations) {
            std::vector<fwdpp::uint_t> rv;
            unsigned nmuts = gsl_ran_poisson(rng.get(), o.mu);
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
    fwdpp::ts::table_index_t first_parental_index = 0,
                             next_index = 2 * pop.diploids.size();
    bool simplified = false;
    std::vector<std::size_t> individual_labels(o.N);
    std::iota(individual_labels.begin(), individual_labels.end(), 0);
    std::vector<std::size_t> individuals;
    if (o.ancient_sample_size > 0)
        {
            individuals.resize(o.ancient_sample_size);
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
    auto ff = fwdpp::multiplicative_diploid(fwdpp::fitness(o.scaling));

    auto genetics = fwdpp::make_genetic_parameters(std::move(ff), std::move(mmodel),
                                                   std::move(recmap));
    auto lookup = calculate_fitnesses(pop, fitnesses, genetics.gvalue);
    std::vector<fwdpp::ts::table_index_t> preserved_nodes;
    for (; generation <= 10 * o.N; ++generation)
        {
            auto pick1 = [&lookup, &rng]() {
                return gsl_ran_discrete(rng.get(), lookup.get());
            };
            auto pick2 = [&lookup, &rng](const std::size_t /*p1*/) {
                return gsl_ran_discrete(rng.get(), lookup.get());
            };
            evolve_generation(rng, pop, genetics, o.N, pick1, pick2, update_offspring,
                              generation, tables, first_parental_index, next_index);
            // Recalculate fitnesses and the lookup table.
            lookup = calculate_fitnesses(pop, fitnesses, genetics.gvalue);
            if (generation % o.gcint == 0.0)
                {
                    auto rv = simplify_tables(
                        pop, generation, pop.mcounts_from_preserved_nodes, tables,
                        simplifier, tables.num_nodes() - 2 * o.N, 2 * o.N,
                        preserved_nodes, o.preserve_fixations);
                    if (!o.preserve_fixations)
                        {
                            genetics.mutation_recycling_bin = fwdpp::ts::make_mut_queue(
                                pop.mcounts, pop.mcounts_from_preserved_nodes);
                        }
                    else
                        {
                            genetics.mutation_recycling_bin = fwdpp::ts::make_mut_queue(
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
                    next_index += 2 * o.N;

                    // TODO: update or remove this bit
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
                    //            tables.mutations.begin(),
                    //            tables.mutations.end(),
                    //            [&pop](const fwdpp::ts::mutation_record &mr) {
                    //                return pop.mcounts[mr.key] == 0;
                    //            });
                    //        tables.mutations.erase(
                    //            itr, tables.mutations.end());
                    //    }
                }
            if (o.ancient_sampling_interval > 0
                && generation % o.ancient_sampling_interval == 0.0
                // This last check forbids us recording the
                // final generation as ancient samples.
                && generation < 10 * o.N)
                {
                    gsl_ran_choose(rng.get(), individuals.data(), individuals.size(),
                                   individual_labels.data(), individual_labels.size(),
                                   sizeof(std::size_t));
                    for (auto i : individuals)
                        {
                            auto x
                                = fwdpp::ts::get_parent_ids(first_parental_index, i, 0);
                            assert(x.first >= first_parental_index);
                            assert(x.second >= first_parental_index);
                            assert(x.first < tables.num_nodes());
                            assert(x.second < tables.num_nodes());
                            assert(tables.nodes[x.first].time == generation);
                            assert(tables.nodes[x.second].time == generation);
                            assert(std::find(preserved_nodes.begin(),
                                             preserved_nodes.end(), x.first)
                                   == preserved_nodes.end());
                            assert(std::find(preserved_nodes.begin(),
                                             preserved_nodes.end(), x.second)
                                   == preserved_nodes.end());
                            preserved_nodes.push_back(x.first);
                            preserved_nodes.push_back(x.second);
                            // Record the metadata for our ancient samples
                            // Here, we record each individual's actual fitness.
                            // If we wanted relative fitness, then
                            // we'd have to copy fitnesses into a temp
                            // and normalize it appropriately.
                            ancient_sample_metadata.emplace_back(
                                i, generation, fitnesses[i], x.first, x.second);
                        }
                }
        }
    if (!simplified)
        {
            auto rv = simplify_tables(pop, generation, pop.mcounts_from_preserved_nodes,
                                      tables, simplifier, tables.num_nodes() - 2 * o.N,
                                      2 * o.N, preserved_nodes, o.preserve_fixations);
            if (o.preserve_fixations)
                {
                    std::vector<std::int32_t> samples(2 * o.N);
                    std::iota(samples.begin(), samples.end(), 0);
                    fwdpp::ts::count_mutations(tables, pop.mutations, samples,
                                               preserved_nodes, pop.mcounts,
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
    // This is an infinite-sites simulation, so we can do some
    // simple tests on our site table:
    if (tables.sites.size() != tables.mutations.size())
        {
            throw std::runtime_error("mutation and site table sizes differ");
        }
    if (!std::is_sorted(begin(tables.sites), end(tables.sites)))
        {
            throw std::runtime_error("site table is not sorted");
        }
    for (std::size_t i = 0; i < tables.sites.size(); ++i)
        {
            if (tables.sites[i].position != pop.mutations[tables.mutations[i].key].pos)
                {
                    throw std::runtime_error("site table data are incorrect");
                }
        }

    // If we have done things correctly, then our
    // ancient sample metadata must match up with
    // what is in our node table.
    for (auto &mr : ancient_sample_metadata)
        {
            if (tables.nodes[mr.n1].time != mr.time
                || tables.nodes[mr.n2].time != mr.time)
                {
                    throw std::runtime_error("invalid ancient sample metadata");
                }
        }
    assert(tables.input_left.size() == tables.edges.size());
    assert(tables.output_right.size() == tables.edges.size());
    std::vector<fwdpp::ts::table_index_t> s(2 * o.N);
    std::iota(s.begin(), s.end(), 0);

    auto neutral_muts
        = apply_neutral_mutations(o, rng, tables, pop, genetics.mutation_recycling_bin);
    if (!std::is_sorted(begin(tables.sites), end(tables.sites)))
        {
            throw std::runtime_error(
                "site table is not sorted after adding neutral mutations");
        }
    for (auto &mr : tables.mutations)
        {
            if (pop.mutations[mr.key].pos != tables.sites[mr.site].position)
                {
                    throw std::runtime_error("site table data incorrect after "
                                             "adding neutral mutations");
                }
        }

    fwdpp::ts::count_mutations(tables, pop.mutations, s, preserved_nodes, pop.mcounts,
                               pop.mcounts_from_preserved_nodes);
    for (std::size_t i = 0; i < pop.mutations.size(); ++i)
        {
            if (pop.mutations[i].neutral)
                {
                    if (!pop.mcounts[i] && !pop.mcounts_from_preserved_nodes[i])
                        {
                            throw std::runtime_error("invalid final mutation count");
                        }
                }
        }
    std::cout << neutral_muts << '\n';

    execute_expensive_leaf_test(o, tables, s, preserved_nodes);
    execute_matrix_test(o, pop, tables, s, preserved_nodes);
    execute_serialization_test(o, tables);
    visit_sites_test(o, pop, tables, s, preserved_nodes);
    write_sfs(o, rng, tables, s);
}
