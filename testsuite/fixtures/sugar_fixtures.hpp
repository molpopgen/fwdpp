#ifndef FWDPP_TESTSUITE_SUGAR_FIXTURES_HPP
#define FWDPP_TESTSUITE_SUGAR_FIXTURES_HPP

#include <fwdpp/popgenmut.hpp>
#include <fwdpp/slocuspop.hpp>
#include <fwdpp/mlocuspop.hpp>
#include <fwdpp/popgenmut.hpp>
#include <fwdpp/GSLrng_t.hpp>
#include <fwdpp/fitness_models.hpp>
#include <fwdpp/extensions/regions.hpp>
#include <fwdpp/simfunctions/recycling.hpp>
#include <testsuite/util/custom_dip.hpp>
#include <functional>
#include <numeric>

struct slocuspop_popgenmut_fixture
{
    using poptype = fwdpp::slocuspop<fwdpp::popgenmut>;
    using rng_t = fwdpp::GSLrng_t<fwdpp::GSL_RNG_TAUS2>;
    poptype pop;
    rng_t rng;
    unsigned generation;
    slocuspop_popgenmut_fixture(const unsigned seed = 0)
        : pop(poptype(1000)), rng(rng_t(seed)), generation(0)
    {
    }
};

struct slocuspop_popgenmut_custom_fixture
{
    using poptype
        = fwdpp::slocuspop<fwdpp::popgenmut, custom_diploid_testing_t>;
    using rng_t = fwdpp::GSLrng_t<fwdpp::GSL_RNG_TAUS2>;
    poptype pop;
    rng_t rng;
    unsigned generation;
    slocuspop_popgenmut_custom_fixture(const unsigned seed = 0)
        : pop(poptype(1000)), rng(rng_t(seed)), generation(0)
    {
    }
};

class mlocuspop_popgenmut_fixture
{
  public:
    using poptype = fwdpp::mlocuspop<fwdpp::popgenmut>;
    using mutmodel = std::function<std::size_t(fwdpp::flagged_mutation_queue &,
                                               poptype::mcont_t &)>;

    static const int nloci;
    static const double per_locus_rate;

  public:
    using rng_t = fwdpp::GSLrng_t<fwdpp::GSL_RNG_TAUS2>;
    using recmodel = std::function<std::vector<double>()>;
    // Fitness function
    struct multilocus_additive
    {
      public:
        using result_type = double;
        inline double
        operator()(const poptype::dipvector_t::value_type &diploid,
                   const poptype::gcont_t &gametes,
                   const poptype::mcont_t &mutations) const
        {
            using dip_t = poptype::dipvector_t::value_type::value_type;
            return std::max(
                0., 1.
                        + std::accumulate(
                              diploid.begin(), diploid.end(), 0.,
                              [&gametes, &mutations](const double d,
                                                     const dip_t &dip) {
                                  return d
                                         + fwdpp::additive_diploid(
                                               fwdpp::fitness(1.))(
                                               gametes[dip.first],
                                               gametes[dip.second], mutations)
                                         - 1.;
                              }));
        }
    };
    poptype pop;
    unsigned generation;
    rng_t rng;
    std::vector<double> mu, rbw;
    std::vector<mutmodel> mutmodels;
    std::vector<recmodel> recmodels;
    std::vector<fwdpp::extensions::discrete_mut_model<poptype::mcont_t>> vdmm;
    std::vector<mutmodel> bound_mmodels;
    std::vector<fwdpp::extensions::discrete_rec_model> vdrm;
    std::vector<std::function<std::vector<double>(void)>> bound_recmodels;
    mlocuspop_popgenmut_fixture(const unsigned seed = 0)
        /*! N=1000, 4 loci */
        : pop(poptype(1000, make_boundaries())), generation(0), rng{ seed },
          mu(std::vector<double>(nloci, per_locus_rate)),
          rbw(std::vector<double>(nloci - 1, per_locus_rate)),
          mutmodels{ make_mutmodels() }, recmodels{ make_recmodels() },
          vdmm(fill_vdmm()),
          bound_mmodels(fwdpp::extensions::bind_vec_dmm(rng.get(), vdmm)),
          vdrm(fill_vdrm()), bound_recmodels{ fill_bound_recmodels() }
    {
    }

  private:
    std::vector<std::function<std::vector<double>(void)>>
    fill_bound_recmodels();
    std::vector<fwdpp::extensions::discrete_mut_model<poptype::mcont_t>>
    fill_vdmm();

    /* We are going to generate a set of recombination
     * regions for a multi-locus
     * simulation.  There will be four loci total.  Each locus
     * (i=0 through 3) will have recombination occurring on
     * the continuous inerval [i,i+1).  Further,
     * each locus will have three regions of different
     * relative recombination rates.  The positions of each
     * region in each locus will be:
     * [i,i+1/3)
     * [i+1/3,i+2/3)
     * [i,i+1)
     * The relative weight on each region will be 1,10,1.
     * The total recombination rate on each region will be 1e-4
     * per diploid, per generation.
     */
    std::vector<fwdpp::extensions::discrete_rec_model> fill_vdrm();

    std::vector<std::pair<double, double>> make_boundaries();

    std::vector<mutmodel> make_mutmodels();
    std::vector<recmodel> make_recmodels();
};

struct slocuspop_objects
{
    using poptype = slocuspop_popgenmut_fixture::poptype;

    poptype::dipvector_t diploids;
    poptype::mcont_t mutations;
    poptype::gcont_t gametes;

    using mutation_container
        = poptype::gcont_t::value_type::mutation_container;

    slocuspop_objects() : diploids{}, mutations{}, gametes{}
    {
        // Add some mutations
        for (unsigned i = 0; i < 3; ++i)
            {
                // position = i, effect size = i
                // meaning muts 1 and 2 not neutral
                mutations.emplace_back(i, i, 1, 0);
            }

        // Add two gametes
        gametes.emplace_back(1, mutation_container{ 0 },
                             mutation_container{ 1 });
        gametes.emplace_back(1, mutation_container{}, mutation_container{ 2 });

        // Add a diploid
        diploids.emplace_back(0, 1);
        assert(gametes.size() == 2);
        assert(mutations.size() == 3);
    }
};

struct mlocuspop_objects
{
    using poptype = mlocuspop_popgenmut_fixture::poptype;

    poptype::dipvector_t diploids;
    poptype::mcont_t mutations;
    poptype::gcont_t gametes;
    decltype(poptype::locus_boundaries) locus_boundaries;

    using mutation_container
        = poptype::gcont_t::value_type::mutation_container;

    mlocuspop_objects()
        : diploids{}, mutations{}, gametes{}, locus_boundaries{}
    {
        // Add some mutations
        for (unsigned i = 0; i < 3; ++i)
            {
                // position = i, effect size = i
                // meaning muts 1 and 2 not neutral
                mutations.emplace_back(i, i, 1, 0);
            }

        // Add two gametes
        gametes.emplace_back(1, mutation_container{ 0 },
                             mutation_container{ 1 });
        gametes.emplace_back(3, mutation_container{}, mutation_container{ 2 });

        diploids.resize(2); // two demes
        // Add a diploid
        diploids[0].emplace_back(0, 1);
        diploids[1].emplace_back(1, 1);
        locus_boundaries.emplace_back(0, 1);
        locus_boundaries.emplace_back(1, 3);
        assert(gametes.size() == 2);
        assert(mutations.size() == 3);
    }
};

#endif
