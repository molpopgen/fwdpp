#ifndef FWDPP_TESTSUITE_SUGAR_FIXTURES_HPP
#define FWDPP_TESTSUITE_SUGAR_FIXTURES_HPP

#include <fwdpp/sugar/popgenmut.hpp>
#include <fwdpp/sugar/slocuspop.hpp>
#include <fwdpp/sugar/mlocuspop.hpp>
#include <fwdpp/sugar/popgenmut.hpp>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/poisson_xover.hpp>
#include <fwdpp/recbinder.hpp>
#include <fwdpp/fitness_models.hpp>
#include <fwdpp/extensions/regions.hpp>
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
    using mutmodel = std::function<std::size_t(std::queue<std::size_t> &,
                                               poptype::mcont_t &)>;

  private:
    std::vector<fwdpp::extensions::discrete_mut_model<poptype::mcont_t>>
    fill_vdmm(std::vector<mutmodel> mutmodels)
    {
        // create a vector of extensions::discrete_mut_model
        std::vector<fwdpp::extensions::discrete_mut_model<poptype::mcont_t>>
            vdmm_;
        for (auto &m : mutmodels)
            {
                vdmm_.emplace_back(
                    fwdpp::extensions::discrete_mut_model<poptype::mcont_t>(
                        std::move(m), std::vector<double>(1, 1.0)));
            }

        return vdmm_;
    }
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
    std::vector<fwdpp::extensions::discrete_rec_model>
    fill_vdrm(const gsl_rng *r)
    {
        double recrate_region = 1e-4;
        // set up a vector of extensions::discrete_rec_model
        // with different regions sizes and weights
        std::vector<fwdpp::extensions::discrete_rec_model> vdrm_;
        for (unsigned i = 0; i < 4; ++i)
            {
                std::vector<fwdpp::extensions::discrete_rec_model::
                                function_type>
                    f;
                std::vector<double> w{ 1., 10., 1. };
                f.push_back([&r, i](std::vector<double> &b) {
                    b.push_back(
                        gsl_ran_flat(r, static_cast<double>(i),
                                     static_cast<double>(i) + 1. / 3.));
                });
                f.push_back([&r, i](std::vector<double> &b) {
                    b.push_back(
                        gsl_ran_flat(r, static_cast<double>(i) + 1. / 3,
                                     static_cast<double>(i) + 2. / 3.));
                });
                f.push_back([&r, i](std::vector<double> &b) {
                    b.push_back(gsl_ran_flat(r, static_cast<double>(i),
                                             static_cast<double>(i) + 1.));
                });
                fwdpp::extensions::discrete_rec_model drm(
                    recrate_region, std::move(f), std::move(w));
                vdrm_.emplace_back(std::move(drm));
            }
        return vdrm_;
    }

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
                0.,
                1. + std::accumulate(diploid.begin(), diploid.end(), 0.,
                                     [&gametes, &mutations](const double d,
                                                            const dip_t &dip) {
                                         return d + fwdpp::additive_diploid()(
                                                        gametes[dip.first],
                                                        gametes[dip.second],
                                                        mutations)
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
        : pop(poptype(1000, 4)), generation(0), rng{ seed },
          mu(std::vector<double>(4, 0.005)),
          rbw(std::vector<double>(3, 0.005)),
          mutmodels(
              { // Locus 0: positions Uniform [0,1)
                [this](std::queue<std::size_t> &recbin,
                       poptype::mcont_t &mutations) {
                    return fwdpp::infsites_popgenmut(
                        recbin, mutations, rng.get(), pop.mut_lookup,
                        this->generation, 0.5,
                        [this]() { return gsl_rng_uniform(rng.get()); },
                        []() { return -0.01; }, []() { return 1.; });
                },
                // Locus 1: positions Uniform [1,2)
                [this](std::queue<std::size_t> &recbin,
                       poptype::mcont_t &mutations) {
                    return fwdpp::infsites_popgenmut(
                        recbin, mutations, rng.get(), pop.mut_lookup,
                        this->generation, 0.5,
                        [this]() { return gsl_ran_flat(rng.get(), 1., 2.); },
                        []() { return -0.01; }, []() { return 1.; });
                },
                // Locus 2: positions Uniform [2,3)
                [this](std::queue<std::size_t> &recbin,
                       poptype::mcont_t &mutations) {
                    return fwdpp::infsites_popgenmut(
                        recbin, mutations, rng.get(), pop.mut_lookup,
                        this->generation, 0.5,
                        [this]() { return gsl_ran_flat(rng.get(), 2., 3.); },
                        []() { return -0.01; }, []() { return 1.; });
                },
                // Locus 3: positions Uniform [3,4)
                [this](std::queue<std::size_t> &recbin,
                       poptype::mcont_t &mutations) {
                    return fwdpp::infsites_popgenmut(
                        recbin, mutations, rng.get(), pop.mut_lookup,
                        this->generation, 0.5,
                        [this]() { return gsl_ran_flat(rng.get(), 1., 2.); },
                        []() { return -0.01; }, []() { return 1.; });
                } }),
          recmodels{ fwdpp::recbinder(fwdpp::poisson_xover(0.005, 0., 1.),
                                      this->rng.get()),
                     fwdpp::recbinder(fwdpp::poisson_xover(0.005, 1., 2.),
                                      this->rng.get()),
                     fwdpp::recbinder(fwdpp::poisson_xover(0.005, 2., 3.),
                                      this->rng.get()),
                     fwdpp::recbinder(fwdpp::poisson_xover(0.005, 3., 4.),
                                      this->rng.get()) },
          vdmm(this->fill_vdmm(mutmodels)),
          bound_mmodels(fwdpp::extensions::bind_vec_dmm(rng.get(), vdmm)),
          vdrm(this->fill_vdrm(rng.get())), bound_recmodels{}
    {
        for (unsigned i = 0; i < vdrm.size(); ++i)
            {
                pop.locus_boundaries.emplace_back(i, i + 1);
                bound_recmodels.emplace_back(
                    [this, i]() { return vdrm.at(i)(rng.get()); });
            }
    }
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

    using mutation_container
        = poptype::gcont_t::value_type::mutation_container;

    mlocuspop_objects() : diploids{}, mutations{}, gametes{}
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
        assert(gametes.size() == 2);
        assert(mutations.size() == 3);
    }
};

#endif
