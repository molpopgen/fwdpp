#ifndef FWDPP_TESTSUITE_SUGAR_FIXTURES_HPP
#define FWDPP_TESTSUITE_SUGAR_FIXTURES_HPP

#include <fwdpp/sugar/popgenmut.hpp>
#include <fwdpp/sugar/singlepop.hpp>
#include <fwdpp/sugar/metapop.hpp>
#include <fwdpp/sugar/multiloc.hpp>
#include <fwdpp/sugar/infsites.hpp>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/recombination.hpp>
#include <fwdpp/fitness_models.hpp>
#include <fwdpp/extensions/regions.hpp>
#include <fwdpp/sugar/serialization.hpp>
#include <testsuite/util/custom_dip.hpp>
#include <functional>

struct singlepop_popgenmut_fixture
{
    using poptype = KTfwd::singlepop<KTfwd::popgenmut>;
    using rng_t = KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2>;
    using mwriter = KTfwd::mutation_writer;
    using mreader = KTfwd::mutation_reader<KTfwd::popgenmut>;
    poptype pop;
    rng_t rng;
    unsigned generation;
    singlepop_popgenmut_fixture(const unsigned seed = 0)
        : pop(poptype(1000)), rng(rng_t(seed)),generation(0)
    {
    }
};

struct singlepop_popgenmut_custom_fixture
{
    using poptype
        = KTfwd::singlepop<KTfwd::popgenmut, custom_diploid_testing_t>;
    using rng_t = KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2>;
    using mwriter = KTfwd::mutation_writer;
    using mreader = KTfwd::mutation_reader<KTfwd::popgenmut>;
    poptype pop;
    rng_t rng;
    unsigned generation;
    singlepop_popgenmut_custom_fixture(const unsigned seed = 0)
        : pop(poptype(1000)), rng(rng_t(seed)),generation(0)
    {
    }
};

struct metapop_popgenmut_fixture
{
    using poptype = KTfwd::metapop<KTfwd::popgenmut>;
    using rng_t = KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2>;
    using mwriter = KTfwd::mutation_writer;
    using mreader = KTfwd::mutation_reader<KTfwd::popgenmut>;
    poptype pop;
    rng_t rng;
    unsigned generation;
    metapop_popgenmut_fixture(const unsigned seed = 0)
        : pop(poptype{ 1000, 1000 }), rng(rng_t(seed)),generation(0)
    {
    }
};

struct metapop_popgenmut_custom_fixture
{
    using poptype = KTfwd::metapop<KTfwd::popgenmut, custom_diploid_testing_t>;
    using rng_t = KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2>;
    using mwriter = KTfwd::mutation_writer;
    using mreader = KTfwd::mutation_reader<KTfwd::popgenmut>;
    poptype pop;
    rng_t rng;
    unsigned generation;
    metapop_popgenmut_custom_fixture(const unsigned seed = 0)
        : pop(poptype{ 1000, 1000 }), rng(rng_t(seed)),generation(0)
    {
    }
};

class multiloc_popgenmut_fixture
{
  private:
    std::vector<KTfwd::extensions::discrete_mut_model>
    fill_vdmm()
    {
        double length = 10.;
        // create a vector of extensions::discrete_mut_model
        std::vector<KTfwd::extensions::discrete_mut_model> vdmm_;
        for (unsigned i = 0; i < 4; ++i)
            {
                double begin = static_cast<double>(i) * length;
                KTfwd::extensions::discrete_mut_model dmm(
                    { begin }, { begin + length }, { 1. }, {}, {}, {}, {});
                vdmm_.emplace_back(std::move(dmm));
            }
        return vdmm_;
    }
    /* We are going to generate a set of recombination
     * regions for a multi-locus
     * simulation.  There will be five loci total.  Each locus
     * (i=0 through 4) will have recombination occurring on
     * the continuous inerval [i*10,i*10+10).  Further,
     * each locus will have three regions of different
     * relative recombination rates.  The positions of each
     * region in each locus will be:
     * [i*10,i*10+3)
     * [i*10+3,i*10+7)
     * [i*10+7,i*10+10)
     * The relative weight on each region will be 1,10,1.
     * The total recombination rate on each region will be 1e-4
     * per diploid, per generation.
     */
    std::vector<KTfwd::extensions::discrete_rec_model>
    fill_vdrm()
    {
        // set up a vector of extensions::discrete_rec_model
        // with different regions sizes and weights
        double length = 10.;
        std::vector<KTfwd::extensions::discrete_rec_model> vdrm_;
        for (unsigned i = 0; i < 4; ++i)
            {
                double begin = static_cast<double>(i) * length;
                KTfwd::extensions::discrete_rec_model drm(
                    { begin, begin + 3., begin + 7. },
                    { begin + 3., begin + 7., begin + length },
                    { 1., 10., 1. });
                vdrm_.emplace_back(std::move(drm));
            }
        return vdrm_;
    }

  public:
    using poptype = KTfwd::multiloc<KTfwd::popgenmut>;
    using rng_t = KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2>;
    using mutmodel = std::function<std::size_t(std::queue<std::size_t> &,
                                               poptype::mcont_t &)>;
    using recmodel = std::function<std::vector<double>(
        const poptype::gamete_t &, const poptype::gamete_t &,
        const poptype::mcont_t &)>;
    using mwriter = KTfwd::mutation_writer;
    using mreader = KTfwd::mutation_reader<KTfwd::popgenmut>;
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
                0., 1. + std::accumulate(
                             diploid.begin(), diploid.end(), 0.,
                             [&gametes, &mutations](const double d,
                                                    const dip_t &dip) {
                                 return d + KTfwd::additive_diploid()(
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
    std::vector<KTfwd::extensions::discrete_mut_model> vdmm;
    std::vector<KTfwd::extensions::discrete_rec_model> vdrm;
    multiloc_popgenmut_fixture(const unsigned seed = 0)
        /*! N=1000, 4 loci */
        : pop(poptype(1000, 4)),
          generation(0),
          rng(rng_t(seed)),
          mu(std::vector<double>(4, 0.005)),
          rbw(std::vector<double>(3, 0.005)),
          mutmodels(
              { // Locus 0: positions Uniform [0,1)
                std::bind(KTfwd::infsites(), std::placeholders::_1,
                          std::placeholders::_2, this->rng.get(),
                          std::ref(pop.mut_lookup), &this->generation, 0.0025,
                          0.0025,
                          [this]() { return gsl_rng_uniform(rng.get()); },
                          []() { return -0.01; }, []() { return 1.; }),
                // Locus 1: positions Uniform [1,2)
                std::bind(KTfwd::infsites(), std::placeholders::_1,
                          std::placeholders::_2, this->rng.get(),
                          std::ref(pop.mut_lookup), &this->generation, 0.0025,
                          0.0025,
                          [this]() { return gsl_ran_flat(rng.get(), 1., 2.); },
                          []() { return -0.01; }, []() { return 1.; }),
                // Locus 2: positions Uniform [2,3)
                std::bind(KTfwd::infsites(), std::placeholders::_1,
                          std::placeholders::_2, this->rng.get(),
                          std::ref(pop.mut_lookup), &this->generation, 0.0025,
                          0.0025,
                          [this]() {
                              return gsl_ran_flat(this->rng.get(), 2., 3.);
                          },
                          []() { return -0.01; }, []() { return 1.; }),
                // Locus 3: positions Uniform [3,4)
                std::bind(KTfwd::infsites(), std::placeholders::_1,
                          std::placeholders::_2, this->rng.get(),
                          std::ref(pop.mut_lookup), &this->generation, 0.0025,
                          0.0025,
                          [this]() { return gsl_ran_flat(rng.get(), 3., 4.); },
                          []() { return -0.01; }, []() { return 1.; }) }),
          recmodels({ std::bind(KTfwd::poisson_xover(), rng.get(), 0.005, 0.,
                                1., std::placeholders::_1,
                                std::placeholders::_2, std::placeholders::_3),
                      std::bind(KTfwd::poisson_xover(), rng.get(), 0.005, 1.,
                                2., std::placeholders::_1,
                                std::placeholders::_2, std::placeholders::_3),
                      std::bind(KTfwd::poisson_xover(), rng.get(), 0.005, 2.,
                                3., std::placeholders::_1,
                                std::placeholders::_2, std::placeholders::_3),
                      std::bind(KTfwd::poisson_xover(), rng.get(), 0.005, 3.,
                                4., std::placeholders::_1,
                                std::placeholders::_2,
                                std::placeholders::_3) }),
          vdmm(this->fill_vdmm()),
          vdrm(this->fill_vdrm())
    {
    }
};

#endif
