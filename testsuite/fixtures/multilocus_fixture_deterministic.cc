#include "multilocus_fixture_deterministic.hpp"

#include <iostream>
#include <vector>
#include <fwdpp/sugar/add_mutation.hpp>
#include <fwdpp/fitness_models.hpp>

using poptype = multilocus_fixture_deterministic::poptype;
const std::size_t multilocus_fixture_deterministic::nloci = 4;
const fwdpp::uint_t multilocus_fixture_deterministic::N = 1000;
const fwdpp::uint_t multilocus_fixture_deterministic::new_mutation_generation
    = 1;

//only valid in conjunction w/make_params2 and swapping NEITHER PARENTAL GAMETE
//p1g1 is mutated at 0.1, 0.41, 1.1, 1.41, 2.1, 2.41, 3.1, 3.41
//p1g2 is mutated at 0.21, 0.61, 1.21, 1.61, 2.21, 2.61, 3.21, 3.61
//p2g1 is mutated at 0.05, 0.45, 1.05, 1.45, 2.05, 2.45, 3.05, 3.45
//p2g2 is mutated at 0.15, 0.75, 1.15, 1.75, 2.15, 2.75, 3.15, 3.75
//p1g1 passes on the following intervals: [0,0.5), [1,1.25), [1.75,2.), [2.5,4.)
//p1g2 passes on the following intervals: [0.5,1), [1.25,1.75), [2,.2.5)
const std::vector<double> multilocus_fixture_deterministic::
    expected_transmitted_mutations_mutate_both_parents_gamete_1
    = {
          0.1,  0.41, 1.1, 3.1, 3.41, // from p1g1
          0.61, 1.61, 2.21            // from p1g2
      };
//p2g1 passes on the following intervals: [0,0.5), [1,1.25), [1.75,2.), [2.5,4.)
//p2g2 passes on the following intervals: [0.5,1), [1.25,1.75), [2,.2.5)
const std::vector<double> multilocus_fixture_deterministic::
    expected_transmitted_mutations_mutate_both_parents_gamete_2
    = {
          0.05, 0.45, 1.05, 3.05, 3.45, //from p2g1
          0.75, 2.15                    // from p2g2
      };

double
multilocus_fixture_deterministic::multilocus_multiplicative::
operator()(const poptype::diploid_t &diploid, const poptype::gcont_t &gametes,
           const poptype::mcont_t &mutations) const
{
    double rv = 1.;
    for (auto &&genotype : diploid)
        {
            rv *= fwdpp::multiplicative_diploid(fwdpp::fitness(2.))(
                genotype, gametes, mutations);
        }
    return std::max(0., rv);
}

void
multilocus_fixture_deterministic::mutate_parent()
{
    std::vector<short> clist(1, 0);
    std::vector<std::size_t> indlist(1, 0);
    for (std::size_t i = 0; i < nloci; ++i)
        {
            // Mutation exactly on left edge of locus
            auto k = fwdpp::add_mutation(pop, i, indlist, clist, i, 0, 0, 0);
            tables.mutation_table.emplace_back(
                fwdpp::ts::mutation_record{ 0, k });
        }
}

void
multilocus_fixture_deterministic::mutate_parent2()
{
    std::vector<short> clist(1, 0);
    std::vector<std::size_t> indlist(1, 0);
    for (std::size_t i = 0; i < nloci; ++i)
        {
            // Mutation exactly on left edge of locus
            auto k = fwdpp::add_mutation(pop, i, indlist, clist, i, 0, 0, 0);
            tables.mutation_table.emplace_back(
                fwdpp::ts::mutation_record{ 0, k });
            k = fwdpp::add_mutation(pop, i, indlist, clist, i + 0.51, 0, 0, 0);
            tables.mutation_table.emplace_back(
                fwdpp::ts::mutation_record{ 0, k });
        }
}

void
multilocus_fixture_deterministic::mutate_both_parents()
// More complex mutations
{
    // Apply het mutations to diploid 0.
    fwdpp::ts::TS_NODE_INT parent_node = 0;
    for (std::size_t i = 0; i < nloci; ++i)
        {
            std::vector<short> clist(1, 0);
            std::vector<std::size_t> indlist(1, 0);
            parent_node = 0;
            auto k = fwdpp::add_mutation(pop, i, indlist, clist, i + 0.1, 0, 0,
                                         0);
            tables.mutation_table.emplace_back(
                fwdpp::ts::mutation_record{ parent_node, k });
            k = fwdpp::add_mutation(pop, i, indlist, clist, i + 0.41, 0, 0, 0);
            tables.mutation_table.emplace_back(
                fwdpp::ts::mutation_record{ parent_node, k });
            // Add mutations to second parental gamete,
            // which corresponds to node 1
            clist = { 1 };
            parent_node = 1;
            k = fwdpp::add_mutation(pop, i, indlist, clist, i + 0.21, 0, 0, 0);
            tables.mutation_table.emplace_back(
                fwdpp::ts::mutation_record{ parent_node, k });
            k = fwdpp::add_mutation(pop, i, indlist, clist, i + 0.61, 0, 0, 0);
            tables.mutation_table.emplace_back(
                fwdpp::ts::mutation_record{ parent_node, k });
        }
    // Add het mutations to diploid 1, whose nodes are 2 and 3.
    for (std::size_t i = 0; i < nloci; ++i)
        {
            parent_node = 2;
            std::vector<short> clist(1, 0);
            std::vector<std::size_t> indlist(1, 1);
            auto k = fwdpp::add_mutation(pop, i, indlist, clist, i + 0.05, 0,
                                         0, 0);
            tables.mutation_table.emplace_back(
                fwdpp::ts::mutation_record{ parent_node, k });
            k = fwdpp::add_mutation(pop, i, indlist, clist, i + 0.45, 0, 0, 0);
            tables.mutation_table.emplace_back(
                fwdpp::ts::mutation_record{ parent_node, k });
            // add mutations to the second gamete of this parent
            clist = { 1 };
            parent_node = 3;
            k = fwdpp::add_mutation(pop, i, indlist, clist, i + 0.15, 0, 0, 0);
            tables.mutation_table.emplace_back(
                fwdpp::ts::mutation_record{ parent_node, k });
            k = fwdpp::add_mutation(pop, i, indlist, clist, i + 0.75, 0, 0, 0);
            tables.mutation_table.emplace_back(
                fwdpp::ts::mutation_record{ parent_node, k });
        }
    // confirm and throw exceptions if we've messed up
    for (int i = 0; i < 2; ++i)
        {
            auto &dip = pop.diploids[i];
            for (auto &locus : dip)
                {
                    if (pop.gametes[locus.first].mutations.size() != 2)
                        {
                            throw std::runtime_error(
                                "first gamete contains wrong number of "
                                "mutations");
                        }
                    if (pop.gametes[locus.second].mutations.size() != 2)
                        {
                            throw std::runtime_error(
                                "second gamete contains wrong number of "
                                "mutations");
                        }
                }
        }
}

auto
multilocus_fixture_deterministic::make_params()
    -> decltype(fwdpp::make_genetic_parameters_with_swapper(
        gvalue, mmodels, intralocus_rec, interlocus_rec, do_not_swap))
// NOTE: we use do_not_swap to suppress any initial randomness
// for the test
{
    return fwdpp::make_genetic_parameters_with_swapper(
        gvalue, mmodels, intralocus_rec, interlocus_rec, do_not_swap);
}

auto
multilocus_fixture_deterministic::make_params_swap_second()
    -> decltype(fwdpp::make_genetic_parameters_with_swapper(
        gvalue, mmodels, intralocus_rec, interlocus_rec, swap_second))
{
    return fwdpp::make_genetic_parameters_with_swapper(
        gvalue, mmodels, intralocus_rec, interlocus_rec,
        swap_second_parent_only());
}

auto
multilocus_fixture_deterministic::make_params2()
    -> decltype(fwdpp::make_genetic_parameters_with_swapper(
        gvalue, mmodels, intralocus_rec2, interlocus_rec2, do_not_swap))
{
    return fwdpp::make_genetic_parameters_with_swapper(
        gvalue, mmodels, intralocus_rec2, interlocus_rec2, do_not_swap);
}

void
multilocus_fixture_deterministic::validate_mutations_positions_1(
    const multilocus_fixture_deterministic::poptype::diploid_t &offspring)
{
    std::vector<fwdpp::uint_t> all_gamete_mut_keys;
    for (auto &locus : offspring)
        {
            all_gamete_mut_keys.insert(
                end(all_gamete_mut_keys),
                begin(pop.gametes[locus.first].mutations),
                end(pop.gametes[locus.first].mutations));
        }
    for (auto p : expected_mutation_positions_1)
        {
            auto itr = std::find_if(begin(all_gamete_mut_keys),
                                    end(all_gamete_mut_keys),
                                    [this, p](fwdpp::uint_t k) {
                                        return pop.mutations[k].pos == p;
                                    });
            if (itr == end(all_gamete_mut_keys))
                {
                    throw std::runtime_error("mutation position not found");
                }
        }
}

void
multilocus_fixture_deterministic::validate_mutations_positions_2(
    const multilocus_fixture_deterministic::poptype::diploid_t &offspring)
{
    std::vector<fwdpp::uint_t> all_gamete_mut_keys;
    for (auto &locus : offspring)
        {
            all_gamete_mut_keys.insert(
                end(all_gamete_mut_keys),
                begin(pop.gametes[locus.first].mutations),
                end(pop.gametes[locus.first].mutations));
        }
    for (auto p : expected_mutation_positions_2)
        {
            auto itr = std::find_if(begin(all_gamete_mut_keys),
                                    end(all_gamete_mut_keys),
                                    [this, p](fwdpp::uint_t k) {
                                        return pop.mutations[k].pos == p;
                                    });
            if (itr == end(all_gamete_mut_keys))
                {
                    throw std::runtime_error("mutation position not found");
                }
        }
}

void
multilocus_fixture_deterministic::validate_mutations_positions_1_recparams2(
    const multilocus_fixture_deterministic::poptype::diploid_t &offspring)
{
    std::vector<fwdpp::uint_t> all_gamete_mut_keys;
    for (auto &locus : offspring)
        {
            all_gamete_mut_keys.insert(
                end(all_gamete_mut_keys),
                begin(pop.gametes[locus.first].mutations),
                end(pop.gametes[locus.first].mutations));
        }
    for (auto p : expected_mutation_positions_1_recparams2)
        {
            auto itr = std::find_if(begin(all_gamete_mut_keys),
                                    end(all_gamete_mut_keys),
                                    [this, p](fwdpp::uint_t k) {
                                        return pop.mutations[k].pos == p;
                                    });
            if (itr == end(all_gamete_mut_keys))
                {
                    throw std::runtime_error("mutation position not found");
                }
        }
}
// PRIVATE FUNCTIONS

std::vector<std::pair<double, double>>
multilocus_fixture_deterministic::make_boundaries()
{
    std::vector<std::pair<double, double>> rv;
    for (std::size_t i = 0; i < nloci; ++i)
        {
            rv.emplace_back(i, i + 1);
        }
    return rv;
}

std::vector<std::function<std::vector<fwdpp::uint_t>(
    fwdpp::flagged_mutation_queue &, poptype::mcont_t &)>>
multilocus_fixture_deterministic::make_mmodels()
// Every locus gets 1 mutation.
{
    std::vector<std::function<std::vector<fwdpp::uint_t>(
        fwdpp::flagged_mutation_queue &, poptype::mcont_t &)>>
        rv;
    const auto s = []() { return -0.01; };
    const auto h = []() { return 1.; };
    for (std::size_t i = 0; i < nloci; ++i)
        {
            const auto generate_mutation_position
                = [this, i]() { return gsl_ran_flat(rng.get(), i, i + 1); };
            const auto make_mutation
                = [this, generate_mutation_position, s,
                   h](fwdpp::flagged_mutation_queue &recbin,
                      poptype::mcont_t &mutations) {
                      return fwdpp::infsites_popgenmut(
                          recbin, mutations, rng.get(), pop.mut_lookup,
                          new_mutation_generation,
                          // 1.0 signifies 100% of mutations will be selected
                          1.0, generate_mutation_position, s, h);
                  };
            const auto mmodel
                = [make_mutation](fwdpp::flagged_mutation_queue &recbin,
                                  poptype::mcont_t &mutations) {
                      std::vector<fwdpp::uint_t> keys;
                      unsigned nmuts = 1;
                      for (unsigned m = 0; m < nmuts; ++m)
                          {
                              keys.push_back(make_mutation(recbin, mutations));
                          }
                      std::sort(begin(keys), end(keys),
                                [&mutations](const fwdpp::uint_t a,
                                             const fwdpp::uint_t b) {
                                    return mutations[a].pos < mutations[b].pos;
                                });
                      return keys;
                  };
            rv.emplace_back(std::move(mmodel));
        }
    return rv;
}

std::vector<std::function<std::vector<double>(void)>>
multilocus_fixture_deterministic::make_intralocus_rec()
{
    std::vector<std::function<std::vector<double>(void)>> rv;
    for (std::size_t i = 0; i < nloci; ++i)
        {
            if (i % 2 == 0.) // No xover
                {
                    rv.emplace_back([]() { return std::vector<double>(); });
                }
            else // xover in middle of locus
                {
                    rv.emplace_back([i]() {
                        return std::vector<double>(
                            { i + 0.5, std::numeric_limits<double>::max() });
                    });
                }
        }
    return rv;
}

std::vector<std::function<std::vector<double>(void)>>
multilocus_fixture_deterministic::make_intralocus_rec2()
{
    std::vector<std::function<std::vector<double>(void)>> rv;
    rv.emplace_back([]() {
        return std::vector<double>(
            { 0.5, std::numeric_limits<double>::max() });
    });
    rv.emplace_back([]() {
        return std::vector<double>(
            { 1.25, 1.75, std::numeric_limits<double>::max() });
    });
    rv.emplace_back([]() {
        return std::vector<double>(
            { 2.5, std::numeric_limits<double>::max() });
    });
    rv.emplace_back([]() { return std::vector<double>(); });
    return rv;
}

std::vector<std::function<unsigned(void)>>
multilocus_fixture_deterministic::make_interlocus_rec()
// The number of recombination events between
// loci i and i + 1 equals i.
{
    std::vector<std::function<unsigned(void)>> rv;
    for (std::size_t i = 0; i < nloci - 1; ++i)
        {
            rv.push_back([i]() -> unsigned { return i; });
        }
    return rv;
}

std::vector<std::function<unsigned(void)>>
multilocus_fixture_deterministic::make_interlocus_rec2()
{
    std::vector<std::function<unsigned(void)>> rv;
    rv.push_back([]() -> unsigned { return 1; });
    rv.push_back([]() -> unsigned { return 1; });
    rv.push_back([]() -> unsigned { return 0; });
    return rv;
}

