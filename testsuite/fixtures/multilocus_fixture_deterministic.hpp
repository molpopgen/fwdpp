#include <fwdpp/ts/table_collection.hpp>
#include <fwdpp/simparams.hpp>
#include <fwdpp/popgenmut.hpp>
#include <fwdpp/mlocuspop.hpp>
#include <fwdpp/GSLrng_t.hpp>
#include <fwdpp/sugar/add_mutation.hpp>

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
    poptype pop;
    fwdpp::ts::table_collection tables;
    fwdpp::GSLrng_mt rng;
    std::vector<std::function<std::vector<fwdpp::uint_t>(
        std::queue<std::size_t> &, poptype::mcont_t &)>>
        mmodels;
    std::vector<std::function<std::vector<double>(void)>> intralocus_rec;
    std::vector<std::function<unsigned(void)>> interlocus_rec;
    std::function<int(const gsl_rng *, std::size_t, std::size_t)> do_not_swap;
    swap_second_parent_only swap_second;
    multilocus_multiplicative gvalue;
    multilocus_fixture_deterministic()
        : pop(N, make_boundaries()), tables(2 * N, 0, 0, nloci), rng(42),
          mmodels(make_mmodels()), intralocus_rec(make_intralocus_rec()),
          interlocus_rec(make_interlocus_rec()),
          do_not_swap(
              [](const gsl_rng *, std::size_t, std::size_t) { return 0; }),
          swap_second(), gvalue()
    {
    }

    // NOTE: we use do_not_swap to suppress any initial randomness
    // for the test
    auto
    make_params() -> decltype(fwdpp::make_genetic_parameters_with_swapper(
        std::move(gvalue), std::move(mmodels), std::move(intralocus_rec),
        std::move(interlocus_rec), do_not_swap));

    // NOTE: we use swap_second to suppress any initial randomness
    // for the test
    auto
    make_params_swap_second()
        -> decltype(fwdpp::make_genetic_parameters_with_swapper(
            std::move(gvalue), std::move(mmodels), std::move(intralocus_rec),
            std::move(interlocus_rec), swap_second));

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

  private:
    std::vector<std::pair<double, double>> make_boundaries();

    // Every locus gets 1 mutation.
    std::vector<std::function<std::vector<fwdpp::uint_t>(
        std::queue<std::size_t> &, poptype::mcont_t &)>>
    make_mmodels();

    std::vector<std::function<std::vector<double>(void)>>
    make_intralocus_rec();

    // The number of recombination events between
    // loci i and i + 1 equals i.
    std::vector<std::function<unsigned(void)>> make_interlocus_rec();
};

