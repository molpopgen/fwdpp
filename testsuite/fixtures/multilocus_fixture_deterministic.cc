#include "multilocus_fixture_deterministic.hpp"

#include <vector>
#include <fwdpp/sugar/add_mutation.hpp>
#include <fwdpp/fitness_models.hpp>

using poptype = multilocus_fixture_deterministic::poptype;
const std::size_t multilocus_fixture_deterministic::nloci = 4;
const fwdpp::uint_t multilocus_fixture_deterministic::N = 1000;

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
auto
multilocus_fixture_deterministic::make_params()
    -> decltype(fwdpp::make_genetic_parameters_with_swapper(
        std::move(gvalue), std::move(mmodels), std::move(intralocus_rec),
        std::move(interlocus_rec), do_not_swap))
// NOTE: we use do_not_swap to suppress any initial randomness
// for the test
{
    return fwdpp::make_genetic_parameters_with_swapper(
        std::move(gvalue), std::move(mmodels), std::move(intralocus_rec),
        std::move(interlocus_rec), do_not_swap);
}

auto
multilocus_fixture_deterministic::make_params_swap_second()
    -> decltype(fwdpp::make_genetic_parameters_with_swapper(
        std::move(gvalue), std::move(mmodels), std::move(intralocus_rec),
        std::move(interlocus_rec), swap_second))
{
    return fwdpp::make_genetic_parameters_with_swapper(
        std::move(gvalue), std::move(mmodels), std::move(intralocus_rec),
        std::move(interlocus_rec), swap_second);
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

std::vector<std::function<std::vector<fwdpp::uint_t>(std::queue<std::size_t> &,
                                                     poptype::mcont_t &)>>
multilocus_fixture_deterministic::make_mmodels()
// Every locus gets 1 mutation.
{
    std::vector<std::function<std::vector<fwdpp::uint_t>(
        std::queue<std::size_t> &, poptype::mcont_t &)>>
        rv;
    const auto s = []() { return -0.01; };
    const auto h = []() { return 1.; };
    for (std::size_t i = 0; i < nloci; ++i)
        {
            const auto generate_mutation_position
                = [this, i]() { return gsl_ran_flat(rng.get(), i, i + 1); };
            const auto make_mutation = [this, generate_mutation_position, s,
                                        h](std::queue<std::size_t> &recbin,
                                           poptype::mcont_t &mutations) {
                return fwdpp::infsites_popgenmut(
                    recbin, mutations, rng.get(), pop.mut_lookup, 0,
                    // 1.0 signifies 100% of mutations will be selected
                    1.0, generate_mutation_position, s, h);
            };
            const auto mmodel
                = [make_mutation](std::queue<std::size_t> &recbin,
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

