#include <config.h>
#include <cmath>
#include <queue>
#include <cstdint>
#include <fwdpp/ts/table_collection.hpp>
#include <fwdpp/ts/generate_offspring.hpp>
#include <fwdpp/simparams.hpp>
#include <fwdpp/popgenmut.hpp>
#include <fwdpp/mlocuspop.hpp>
#include <fwdpp/fitness_models.hpp>
#include <fwdpp/GSLrng_t.hpp>
#include <boost/test/unit_test.hpp>

struct multilocus_fixture_deterministic
// Deterministically adds mutations and recombination
// events to a multi-locus population
{

  public:
    using poptype = fwdpp::mlocuspop<fwdpp::popgenmut>;
    struct multilocus_multiplicative
    {
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
    multilocus_multiplicative gvalue;
    multilocus_fixture_deterministic()
        : pop(N, make_boundaries()), tables(2 * N, 0, 0, nloci), rng(42),
          mmodels(make_mmodels()), intralocus_rec(make_intralocus_rec()),
          interlocus_rec(make_interlocus_rec()), gvalue()
    {
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

BOOST_FIXTURE_TEST_CASE(check_multilocus_deterministic_fixture,
                        multilocus_fixture_deterministic)
{
    BOOST_REQUIRE_EQUAL(tables.genome_length(), nloci);
    BOOST_REQUIRE_EQUAL(pop.locus_boundaries.size(), nloci);
    BOOST_REQUIRE_EQUAL(mmodels.size(), nloci);
    BOOST_REQUIRE_EQUAL(intralocus_rec.size(), nloci);
    BOOST_REQUIRE_EQUAL(interlocus_rec.size(), nloci - 1);

    auto params = fwdpp::make_genetic_parameters(
        std::move(gvalue), std::move(mmodels), std::move(intralocus_rec),
        std::move(interlocus_rec));
    for (std::size_t i = 0; i < nloci; ++i)
        {
            auto k = params.generate_mutations[i](
                params.mutation_recycling_bin, pop.mutations);
            BOOST_REQUIRE_EQUAL(k.size(), 1);
            BOOST_REQUIRE_EQUAL(pop.mutations[k[0]].pos >= i, true);
            BOOST_REQUIRE_EQUAL(pop.mutations[k[0]].pos < i + 1, true);
            auto r = params.generate_breakpoints[i]();
            if (i % 2 == 0.)
                {
                    BOOST_REQUIRE_EQUAL(r.size(), 0);
                }
            else
                {
                    BOOST_REQUIRE_EQUAL(r.size(), 2);
                    BOOST_REQUIRE_EQUAL(r[0], i + 0.5);
                }
            if (i)
                {
                    auto ir = params.interlocus_recombination[i - 1]();
                    BOOST_REQUIRE_EQUAL(ir, i - 1);
                }
        }
}

BOOST_FIXTURE_TEST_CASE(test_multilocus_determinisic,
                        multilocus_fixture_deterministic)
{
    auto params = fwdpp::make_genetic_parameters(
        std::move(gvalue), std::move(mmodels), std::move(intralocus_rec),
        std::move(interlocus_rec));
    poptype::diploid_t offspring;
    auto data_to_record = fwdpp::ts::generate_offspring(
        rng.get(), std::make_pair(0, 1), fwdpp::ts::selected_variants_only(),
        pop, params, offspring);
    BOOST_REQUIRE_EQUAL(offspring.size(), nloci);
    for (auto &locus : offspring)
        {
            BOOST_CHECK_EQUAL(pop.gametes[locus.first].mutations.size(), 0);
            BOOST_CHECK_EQUAL(pop.gametes[locus.first].smutations.size(), 1);
            BOOST_CHECK_EQUAL(pop.gametes[locus.second].mutations.size(), 0);
            BOOST_CHECK_EQUAL(pop.gametes[locus.second].smutations.size(),
                                1);
        }
    BOOST_CHECK_EQUAL(data_to_record.first.mutation_keys.size(), nloci);
    BOOST_CHECK_EQUAL(data_to_record.second.mutation_keys.size(), nloci);
}
