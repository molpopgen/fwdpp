#include <fwdpp/ts/table_collection.hpp>
#include <fwdpp/simparams.hpp>
#include <fwdpp/popgenmut.hpp>
#include <fwdpp/mlocuspop.hpp>
#include <fwdpp/fitness_models.hpp>
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
        inline double
        operator()(const poptype::diploid_t &diploid,
                   const poptype::gcont_t &gametes,
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

    auto
    make_params() -> decltype(fwdpp::make_genetic_parameters_with_swapper(
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
    make_params_swap_second()
        -> decltype(fwdpp::make_genetic_parameters_with_swapper(
            std::move(gvalue), std::move(mmodels), std::move(intralocus_rec),
            std::move(interlocus_rec), swap_second))
    // NOTE: we use swap_second to suppress any initial randomness
    // for the test
    {
        return fwdpp::make_genetic_parameters_with_swapper(
            std::move(gvalue), std::move(mmodels), std::move(intralocus_rec),
            std::move(interlocus_rec), swap_second);
    }

  private:
    std::vector<std::pair<double, double>>
    make_boundaries()
    {
        std::vector<std::pair<double, double>> rv;
        for (std::size_t i = 0; i < nloci; ++i)
            {
                rv.emplace_back(i, i + 1);
            }
        return rv;
    }

    std::vector<std::function<std::vector<fwdpp::uint_t>(
        std::queue<std::size_t> &, poptype::mcont_t &)>>
    make_mmodels()
    // Every locus gets 1 mutation.
    {
        std::vector<std::function<std::vector<fwdpp::uint_t>(
            std::queue<std::size_t> &, poptype::mcont_t &)>>
            rv;
        const auto s = []() { return -0.01; };
        const auto h = []() { return 1.; };
        for (std::size_t i = 0; i < nloci; ++i)
            {
                const auto generate_mutation_position = [this, i]() {
                    return gsl_ran_flat(rng.get(), i, i + 1);
                };
                const auto make_mutation = [this, generate_mutation_position,
                                            s,
                                            h](std::queue<std::size_t> &recbin,
                                               poptype::mcont_t &mutations) {
                    return fwdpp::infsites_popgenmut(
                        recbin, mutations, rng.get(), pop.mut_lookup, 0,
                        // 1.0 signifies 100% of mutations will be selected
                        1.0, generate_mutation_position, s, h);
                };
                const auto mmodel = [make_mutation](
                                        std::queue<std::size_t> &recbin,
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
    make_intralocus_rec()
    {
        std::vector<std::function<std::vector<double>(void)>> rv;
        for (std::size_t i = 0; i < nloci; ++i)
            {
                if (i % 2 == 0.) // No xover
                    {
                        rv.emplace_back(
                            []() { return std::vector<double>(); });
                    }
                else // xover in middle of locus
                    {
                        rv.emplace_back([i]() {
                            return std::vector<double>(
                                { i + 0.5,
                                  std::numeric_limits<double>::max() });
                        });
                    }
            }
        return rv;
    }

    std::vector<std::function<unsigned(void)>>
    make_interlocus_rec()
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
};

const std::size_t multilocus_fixture_deterministic::nloci = 4;
const fwdpp::uint_t multilocus_fixture_deterministic::N = 1000;
