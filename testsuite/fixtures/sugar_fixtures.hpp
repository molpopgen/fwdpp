#ifndef FWDPP_TESTSUITE_SUGAR_FIXTURES_HPP
#define FWDPP_TESTSUITE_SUGAR_FIXTURES_HPP

#include <fwdpp/popgenmut.hpp>
#include <fwdpp/diploid_population.hpp>
#include <fwdpp/popgenmut.hpp>
#include <fwdpp/GSLrng_t.hpp>
#include <fwdpp/fitness_models.hpp>
#include <fwdpp/extensions/regions.hpp>
#include <fwdpp/simfunctions/recycling.hpp>
#include <testsuite/util/custom_dip.hpp>
#include <functional>
#include <numeric>

struct diploid_population_popgenmut_fixture
{
    using poptype = fwdpp::diploid_population<fwdpp::popgenmut>;
    using rng_t = fwdpp::GSLrng_t<fwdpp::GSL_RNG_TAUS2>;
    poptype pop;
    rng_t rng;
    unsigned generation;
    diploid_population_popgenmut_fixture(const unsigned seed = 0)
        : pop(poptype(1000)), rng(rng_t(seed)), generation(0)
    {
    }
};

struct diploid_population_popgenmut_custom_fixture
{
    using poptype = fwdpp::diploid_population<fwdpp::popgenmut,
                                              custom_diploid_testing_t>;
    using rng_t = fwdpp::GSLrng_t<fwdpp::GSL_RNG_TAUS2>;
    poptype pop;
    rng_t rng;
    unsigned generation;
    diploid_population_popgenmut_custom_fixture(const unsigned seed = 0)
        : pop(poptype(1000)), rng(rng_t(seed)), generation(0)
    {
    }
};

struct diploid_population_objects
{
    using poptype = diploid_population_popgenmut_fixture::poptype;

    poptype::dipvector_t diploids;
    poptype::mcont_t mutations;
    poptype::gcont_t haploid_genomes;

    using mutation_container
        = poptype::gcont_t::value_type::mutation_container;

    diploid_population_objects() : diploids{}, mutations{}, haploid_genomes{}
    {
        // Add some mutations
        for (unsigned i = 0; i < 3; ++i)
            {
                // position = i, effect size = i
                // meaning muts 1 and 2 not neutral
                mutations.emplace_back(i, i, 1, 0);
            }

        // Add two haploid_genomes
        haploid_genomes.emplace_back(1, mutation_container{ 0 },
                                     mutation_container{ 1 });
        haploid_genomes.emplace_back(1, mutation_container{},
                                     mutation_container{ 2 });

        // Add a diploid
        diploids.emplace_back(0, 1);
        assert(haploid_genomes.size() == 2);
        assert(mutations.size() == 3);
    }
};

#endif
