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
#include <fwdpp/ts/recycling.hpp>
#include <fwdpp/ts/mutate_tables.hpp>
#include <fwdpp/ts/generate_data_matrix.hpp>
#include <fwdpp/ts/remove_fixations_from_gametes.hpp>
#include <fwdpp/GSLrng_t.hpp>
#include <fwdpp/popgenmut.hpp>
#include <fwdpp/slocuspop.hpp>
#include <fwdpp/util.hpp>
#include <fwdpp/fitness_models.hpp>
#include <fwdpp/extensions/callbacks.hpp>
#include <fwdpp/poisson_xover.hpp>
#include <fwdpp/recbinder.hpp>
#include <fwdpp/simparams.hpp>

#include "simplify_tables.hpp"
#include "evolve_generation_ts.hpp"
#include "calculate_fitnesses.hpp"
#include "tree_sequence_examples_common.hpp"

namespace po = boost::program_options;
using poptype = fwdpp::slocuspop<fwdpp::popgenmut>;
using GSLrng = fwdpp::GSLrng_mt;

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

    poptype pop(o.N);
    fwdpp::ts::table_collection tables(2 * pop.diploids.size(), 0, 0, 1.0);
    fwdpp::ts::table_simplifier simplifier(1.0);
    unsigned generation = 1;
    double recrate = o.rho / static_cast<double>(4 * o.N);
    const auto recmap
        = fwdpp::recbinder(fwdpp::poisson_xover(recrate, 0., 1.), rng.get());

    auto get_selection_coefficient
        = make_dfe(o.N, rng, o.mean, o.shape, o.scoeff);
    const auto generate_mutation_position
        = [&rng]() { return gsl_rng_uniform(rng.get()); };
    const auto generate_h = [&o]() { return o.dominance; };
    const auto make_mutation
        = [&pop, &rng, &generation, generate_mutation_position,
           get_selection_coefficient, generate_h](
              std::queue<std::size_t> &recbin, poptype::mcont_t &mutations) {
              return fwdpp::infsites_popgenmut(
                  recbin, mutations, rng.get(), pop.mut_lookup, generation,
                  // 1.0 signifies 100% of mutations will be selected
                  1.0, generate_mutation_position, get_selection_coefficient,
                  generate_h);
          };

    const auto mmodel = [&rng, &o,
                         &make_mutation](std::queue<std::size_t> &recbin,
                                         poptype::mcont_t &mutations) {
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
    fwdpp::ts::TS_NODE_INT first_parental_index = 0,
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

    auto genetics = fwdpp::make_genetic_parameters(
        std::move(ff), std::move(mmodel), std::move(recmap));
    auto lookup = calculate_fitnesses(pop, fitnesses, genetics.gvalue);
    for (; generation <= 10 * o.N; ++generation)
        {
            auto pick1 = [&lookup, &rng]() {
                return gsl_ran_discrete(rng.get(), lookup.get());
            };
            auto pick2 = [&lookup, &rng](const std::size_t /*p1*/) {
                return gsl_ran_discrete(rng.get(), lookup.get());
            };
            evolve_generation(rng, pop, genetics, o.N, pick1, pick2,
                              update_offspring, generation, tables,
                              first_parental_index, next_index);
            // Recalculate fitnesses and the lookup table.
            lookup = calculate_fitnesses(pop, fitnesses, genetics.gvalue);
            if (generation % o.gcint == 0.0)
                {
                    auto rv = simplify_tables(
                        pop, generation, pop.mcounts_from_preserved_nodes,
                        tables, simplifier, tables.num_nodes() - 2 * o.N,
                        2 * o.N, o.preserve_fixations);
                    if (!o.preserve_fixations)
                        {
                            genetics.mutation_recycling_bin
                                = fwdpp::ts::make_mut_queue(
                                    pop.mcounts,
                                    pop.mcounts_from_preserved_nodes);
                        }
                    else
                        {
                            genetics.mutation_recycling_bin
                                = fwdpp::ts::make_mut_queue(
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
                    //            tables.mutation_table.begin(),
                    //            tables.mutation_table.end(),
                    //            [&pop](const fwdpp::ts::mutation_record &mr) {
                    //                return pop.mcounts[mr.key] == 0;
                    //            });
                    //        tables.mutation_table.erase(
                    //            itr, tables.mutation_table.end());
                    //    }
                }
            if (o.ancient_sampling_interval > 0
                && generation % o.ancient_sampling_interval == 0.0
                // This last check forbids us recording the
                // final generation as ancient samples.
                && generation < 10 * o.N)
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
                                      simplifier, tables.num_nodes() - 2 * o.N,
                                      2 * o.N, o.preserve_fixations);
            if (o.preserve_fixations)
                {
                    std::vector<std::int32_t> samples(2 * o.N);
                    std::iota(samples.begin(), samples.end(), 0);
                    fwdpp::ts::count_mutations(
                        tables, pop.mutations, samples, pop.mcounts,
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
    std::vector<fwdpp::ts::TS_NODE_INT> s(2 * o.N);
    std::iota(s.begin(), s.end(), 0);
    const auto neutral_variant_maker
        = [&rng, &pop, &genetics](const double left, const double right,
                                  const fwdpp::uint_t generation) {
              return fwdpp::infsites_popgenmut(
                  genetics.mutation_recycling_bin, pop.mutations, rng.get(),
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
    std::cout << neutral_muts << '\n';

    if (o.leaf_test)
        {
            std::cerr << "Starting sample list validation.  This may take a "
                         "while!\n";
            expensive_leaf_test(tables, s);
            std::cout << "Passed with respect to last generation.\n";
            expensive_leaf_test(tables, tables.preserved_nodes);
            std::cout << "Passed with respect to preserved samples.\n";
        }

    if (o.matrix_test)
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
    if (!o.filename.empty())
        {
            test_serialization(tables, o.filename);
        }

    if (!o.sfsfilename.empty())
        {
            if (!(o.nsam > 2))
                {
                    throw std::invalid_argument(
                        "sample size for site frequency spectrum must be > 2");
                }
            // Simplify w.r.to 100 samples
            std::vector<fwdpp::ts::TS_NODE_INT> small_sample(o.nsam);
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
            std::ofstream sfs_stream(o.sfsfilename.c_str());
            for (std::size_t i = 0; i < sfs.size(); ++i)
                {
                    sfs_stream << (i + 1) << ' ' << sfs[i] << '\n';
                }
        }
}
