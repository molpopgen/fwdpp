/*! \file sugar_metapop_custom_diploid.cc
  \ingroup unit
  \brief Testing fwdpp::metapop with custom diploid type
*/
#include <config.h>
#include <memory>
#include <algorithm>
#include <boost/test/unit_test.hpp>
#include <fwdpp/diploid.hh>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/sugar/singlepop.hpp>
#include <fwdpp/sugar/metapop.hpp>
#include <fwdpp/sugar/infsites.hpp>
#include <fwdpp/sugar/serialization.hpp>
#include <testsuite/util/custom_dip.hpp>
#include <testsuite/util/migpop.hpp>

using mutation_t = fwdpp::popgenmut;

using spoptype = fwdpp::singlepop<mutation_t, custom_diploid_testing_t>;
using poptype = fwdpp::metapop<mutation_t, custom_diploid_testing_t>;

/*
  These next two derived classes mimic what software
  packages like fwdpy do, which is to extend sugar types
  with data that they need.
*/
struct spop_derived : public spoptype
{
    unsigned generation;
    spop_derived(unsigned N) : spoptype(N), generation(0) {}
};

struct mpop_derived : public poptype
{
    unsigned generation;
    mpop_derived(std::initializer_list<unsigned> Ns)
        : poptype(Ns), generation(0)
    {
    }
    mpop_derived(const spop_derived &s) : poptype(s), generation(s.generation)
    {
    }
    mpop_derived(spop_derived &&s)
        : poptype(std::move(s)), generation(s.generation)
    {
    }
};

void
simulate(poptype &pop)
{
    // Evolve for 10 generations
    std::vector<std::function<double(const poptype::diploid_t &,
                                     const poptype::gcont_t &,
                                     const poptype::mcont_t &)>>
        fitness_funcs(2, [](const poptype::diploid_t &d,
                            const poptype::gcont_t &g,
                            const poptype::mcont_t &m) {
            return fwdpp::multiplicative_diploid()(g[d.first], g[d.second], m);
        });
    fwdpp::GSLrng_t<fwdpp::GSL_RNG_TAUS2> rng(0u);
    for (unsigned generation = 0; generation < 10; ++generation)
        {
            std::vector<double> wbar = fwdpp::sample_diploid(
                rng.get(), pop.gametes, pop.diploids, pop.mutations,
                pop.mcounts, &pop.Ns[0], 0.005,
                std::bind(fwdpp::infsites(), std::placeholders::_1,
                          std::placeholders::_2, rng.get(),
                          std::ref(pop.mut_lookup), generation, 0.005, 0.,
                          [&rng]() { return gsl_rng_uniform(rng.get()); },
                          []() { return 0.; }, []() { return 0.; }),
                fwdpp::poisson_xover(rng.get(), 0.005, 0., 1.), fitness_funcs,
                std::bind(migpop, std::placeholders::_1, rng.get(), 0.001),
                pop.neutral, pop.selected);
            fwdpp::update_mutations(pop.mutations, pop.fixations,
                                    pop.fixation_times, pop.mut_lookup,
                                    pop.mcounts, generation, 4000);
        }
}

BOOST_AUTO_TEST_CASE(metapop_sugar_custom_test1)
{
    poptype pop({ 1000, 1000 });
    simulate(pop);
    poptype pop2(pop);
    BOOST_CHECK_EQUAL(pop == pop2, true);
}

BOOST_AUTO_TEST_CASE(metapop_sugar_custom_test2)
{
    poptype pop({ 1000, 1000 });
    simulate(pop);
    for (unsigned i = 0; i < pop.diploids[0].size(); ++i)
        {
            pop.diploids[0][i].i = i;
        }
    poptype pop2{ 0, 0 };
    fwdpp::serialize s;
    std::stringstream buffer;
    s(buffer, pop);
    fwdpp::deserialize()(pop2, buffer);
    BOOST_CHECK_EQUAL(pop == pop2, true);
}

BOOST_AUTO_TEST_CASE(metapop_sugar_custom_test3)
{
    poptype pop({ 1000, 1000 });
    simulate(pop);
    poptype pop2(std::move(pop));
    BOOST_CHECK_EQUAL(pop == pop2, false);
}

BOOST_AUTO_TEST_CASE(metapop_sugar_custom_test4)
{
    poptype pop({ 1000, 1000 });
    simulate(pop);
    poptype pop2 = std::move(pop);
    BOOST_CHECK_EQUAL(pop == pop2, false);
}

// Test construction from singlepops
BOOST_AUTO_TEST_CASE(copy_construct_spop_derived)
{
    spop_derived pop(1000);
    BOOST_REQUIRE_EQUAL(pop.N, 1000);
    mpop_derived mpop(pop);
    BOOST_REQUIRE_EQUAL(mpop.Ns.size(), 1);
    BOOST_REQUIRE_EQUAL(mpop.Ns[0], 1000);
    BOOST_REQUIRE_EQUAL(mpop.diploids.size(), 1);
    BOOST_REQUIRE_EQUAL(mpop.diploids[0].size(), 1000);
    BOOST_REQUIRE_EQUAL(pop.mutations == mpop.mutations, true);
    BOOST_REQUIRE_EQUAL(pop.gametes == mpop.gametes, true);
    BOOST_REQUIRE_EQUAL(pop.diploids == mpop.diploids[0], true);
}

BOOST_AUTO_TEST_CASE(move_construct_spop_derived)
{
    spop_derived pop(1000);
    BOOST_REQUIRE_EQUAL(pop.N, 1000);
    mpop_derived mpop(std::move(pop));
    BOOST_REQUIRE_EQUAL(mpop.Ns.size(), 1);
    BOOST_REQUIRE_EQUAL(mpop.Ns[0], 1000);
    BOOST_REQUIRE_EQUAL(mpop.diploids.size(), 1);
    BOOST_REQUIRE_EQUAL(mpop.diploids[0].size(), 1000);
}

// This is what fwdpy does internally: move stuff in and around shared_ptr
// objects
BOOST_AUTO_TEST_CASE(move_construct_spop_derived_shared_ptr)
{
    std::shared_ptr<spop_derived> pop(new spop_derived(1000));
    BOOST_REQUIRE_EQUAL(pop->N, 1000);
    std::shared_ptr<mpop_derived> mpop(
        new mpop_derived(std::initializer_list<unsigned>{ 0 }));
    mpop.reset(new mpop_derived(*pop.get()));
    BOOST_REQUIRE_EQUAL(mpop->Ns.size(), 1);
    BOOST_REQUIRE_EQUAL(mpop->Ns[0], 1000);
    BOOST_REQUIRE_EQUAL(mpop->diploids.size(), 1);
    BOOST_REQUIRE_EQUAL(mpop->diploids[0].size(), 1000);
}
