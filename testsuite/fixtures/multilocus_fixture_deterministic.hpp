#ifndef FWDPP_TESTSUITE_MULTILOCUS_FIXTURES_DETERMINISTIC_HPP
#define FWDPP_TESTSUITE_MULTILOCUS_FIXTURES_DETERMINISTIC_HPP

#include <fwdpp/ts/table_collection.hpp>
#include <fwdpp/simfunctions/recycling.hpp>
#include <fwdpp/simparams.hpp>
#include <fwdpp/popgenmut.hpp>
#include <fwdpp/mlocuspop.hpp>
#include <fwdpp/GSLrng_t.hpp>

struct multilocus_fixture_deterministic
// Deterministically adds mutations and recombination
// events to a multi-locus population
// The number of cross-overs w/in and b/w loci look like this:
// Locus Within Between*
// 0     0      NA
// 1     1      0
// 2     0      1
// 3     1      2
//
// *Between refers to the number of crossovers occuring between
// loci i-1 and i.
//
// When a crossover happens w/in locus i, it happens at position i+0.5.
//
// Every locus gets 1 mutation uniformly on [i,i+1).  This has to be random,
// so that we can use existing infinite-sites mutation code.  We fix the RNG
// seed.
{

  public:
    using poptype = fwdpp::mlocuspop<fwdpp::popgenmut>;
    struct multilocus_multiplicative
    {
        //NOTE: older compilers may not accept initialization
        //of this type into fwdpp::simparams unless this default
        //constructor is defined:
        multilocus_multiplicative() {}
        double operator()(const poptype::diploid_t &diploid,
                          const poptype::gcont_t &gametes,
                          const poptype::mcont_t &mutations) const;
    };

    class swap_second_parent_only
    {
      private:
        mutable int ncalls;

      public:
        swap_second_parent_only() : ncalls(0) {}
        inline int
        operator()(const gsl_rng *, std::size_t, std::size_t) const
        {
            int rv = ncalls;
            ncalls++;
            return rv;
        }
    };

    static const std::size_t nloci;
    static const fwdpp::uint_t N;
    static const fwdpp::uint_t new_mutation_generation;
    poptype pop;
    fwdpp::ts::table_collection tables;
    fwdpp::GSLrng_mt rng;
    std::vector<std::function<std::vector<fwdpp::uint_t>(
        fwdpp::flagged_mutation_queue &, poptype::mcont_t &)>>
        mmodels;
    std::vector<std::function<std::vector<double>(void)>> intralocus_rec,
        intralocus_rec2;
    std::vector<std::function<unsigned(void)>> interlocus_rec, interlocus_rec2;
    std::function<int(const gsl_rng *, std::size_t, std::size_t)> do_not_swap;
    swap_second_parent_only swap_second;
    multilocus_multiplicative gvalue;
    std::vector<double> expected_breakpoints;
    std::vector<double> expected_breakpoints2;
    // The expected mutation positions are what we expect
    // after calls to mutate_parent and mutate_parent2,
    // respectively:
    std::vector<double> expected_mutation_positions_1,
        expected_mutation_positions_2,
        expected_mutation_positions_1_recparams2;
    static const std::vector<double>
        expected_transmitted_mutations_mutate_both_parents_gamete_1,
        expected_transmitted_mutations_mutate_both_parents_gamete_2;
    decltype(fwdpp::make_genetic_parameters_with_swapper(
        gvalue, mmodels, intralocus_rec, interlocus_rec,
        do_not_swap)) params_no_swap;
    decltype(fwdpp::make_genetic_parameters_with_swapper(
        gvalue, mmodels, intralocus_rec, interlocus_rec,
        swap_second)) params_swap_second;
    decltype(fwdpp::make_genetic_parameters_with_swapper(
        gvalue, mmodels, intralocus_rec2, interlocus_rec2,
        do_not_swap)) params_no_swap2;

    multilocus_fixture_deterministic()
        : pop(N, make_boundaries()), tables(2 * N, 0, 0, nloci), rng(42),
          mmodels(make_mmodels()), intralocus_rec(make_intralocus_rec()),
          intralocus_rec2(make_intralocus_rec2()),
          interlocus_rec(make_interlocus_rec()),
          interlocus_rec2(make_interlocus_rec2()),
          do_not_swap(
              [](const gsl_rng *, std::size_t, std::size_t) { return 0; }),
          swap_second(), gvalue(),
          expected_breakpoints{ { 1.5, 2., 3.5,
                                  std::numeric_limits<double>::max() } },
          expected_breakpoints2({ 0.5, 1., 1.25, 1.75, 2., 2.5,
                                  std::numeric_limits<double>::max() }),
          expected_mutation_positions_1{ { 1., 2., 3., 0. } },
          expected_mutation_positions_2{ { 1., 2., 3., 0., 0.51, 2.51 } },
          expected_mutation_positions_1_recparams2{ { 0., 1., 2.51, 3.,
                                                      3.51 } },
          params_no_swap(make_params()),
          params_swap_second(make_params_swap_second()),
          params_no_swap2(make_params2())
    {
    }

    // We add a variant to each gamete that
    // is exactly at the start of each locus
    // This happens on diploid 0's first gamete
    // at all loci
    void mutate_parent();

    // We add a variant to each gamete that
    // is exactly at the start of each locus,
    // and one mutation within each locus
    // This happens on diploid 0's first gamete
    // at all loci.
    void mutate_parent2();

    // Mutate both chromosomes of both parents
    void mutate_both_parents();

    void validate_mutations_positions_1(const poptype::diploid_t &offspring);

    void validate_mutations_positions_2(const poptype::diploid_t &offspring);

    void validate_mutations_positions_1_recparams2(
        const poptype::diploid_t &offspring);

  private:
    // NOTE: we use do_not_swap to suppress any initial randomness
    // for the test
    auto make_params() -> decltype(fwdpp::make_genetic_parameters_with_swapper(
        gvalue, mmodels, intralocus_rec, interlocus_rec, do_not_swap));

    // NOTE: we use swap_second to suppress any initial randomness
    // for the test
    auto make_params_swap_second()
        -> decltype(fwdpp::make_genetic_parameters_with_swapper(
            gvalue, mmodels, intralocus_rec, interlocus_rec, swap_second));

    auto make_params2()
        -> decltype(fwdpp::make_genetic_parameters_with_swapper(
            gvalue, mmodels, intralocus_rec2, interlocus_rec2, do_not_swap));
    std::vector<std::pair<double, double>> make_boundaries();

    // Every locus gets 1 mutation.  These occur at random
    // positions within each locus
    std::vector<std::function<std::vector<fwdpp::uint_t>(
        fwdpp::flagged_mutation_queue &, poptype::mcont_t &)>>
    make_mmodels();

    std::vector<std::function<std::vector<double>(void)>>
    make_intralocus_rec();

    std::vector<std::function<std::vector<double>(void)>>
    make_intralocus_rec2();

    // The number of recombination events between
    // loci i and i + 1 equals i.
    std::vector<std::function<unsigned(void)>> make_interlocus_rec();
    std::vector<std::function<unsigned(void)>> make_interlocus_rec2();
};

#endif
