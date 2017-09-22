/*
  Separate sex model with sex-specific mutation rates
  and fitness effects of mutations being different in each sex.

  \include sex_limited.cc
  Total madness:
  1. Different mutation rates to variants affecting traits in each sex
  2. Mutation effects are sex-limited ("male" mutations only affect trait
  values in males, etc.)
  3. House of cards model of additive fitness effects and then stabilizing
  selection on trait according to a unit Gaussian
*/

#include <config.h>

#include <limits>
#include <algorithm>
#include <iostream>

// Main fwdpp library header
#include <fwdpp/diploid.hh>
#ifdef HAVE_LIBSEQUENCE
#include <Sequence/SimData.hpp>
#endif
// Include the necessary "sugar" components
#include <fwdpp/sugar/singlepop.hpp>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/sugar/serialization.hpp>
#include <fwdpp/sugar/infsites.hpp>
#include <fwdpp/experimental/sample_diploid.hpp>
// FWDPP-related stuff

struct sex_specific_mutation : public KTfwd::mutation_base
{
    // Effect size
    double s;
    // The effect size applies to this sex, is 0 otherwise
    bool sex;

    sex_specific_mutation(const double &__pos, const double &__s,
                          const bool &__sex, const bool &__neutral)
        : KTfwd::mutation_base(__pos, __neutral), s(__s), sex(__sex)
    {
    }
};

using mtype = sex_specific_mutation;

// We need to define a custom diploid genotype for our model
struct diploid_t
{
    using first_type = std::size_t;
    using second_type = std::size_t;
    first_type first;
    second_type second;
    bool sex; // male = false, I guess...
    // constructors, etc.
    diploid_t() : first(first_type()), second(second_type()), sex(false) {}
    //"perfect forwarding" constructor does not work with iterator from boost
    // containers...
    // diploid_t(first_type && g1, first_type && g2) :
    // first(std::forward(g1)),second(std::forward(g2)),i(numeric_limits<unsigned>::max())
    // {}
    diploid_t(first_type g1, first_type g2) : first(g1), second(g2), sex(false)
    {
    }
    // The following constructors SHOULD be generated automagically by your
    // compiler, so you don't have to:
    //(no idea what, if any, performance effect this may have.  Worst case is
    // prob. the move constructor doesn't get auto-generated...
    // diploid_t( const diploid_t & ) = default;
    // diploid_t( diploid_t && ) = default;
    // diploid_t & operator=(const diploid_t &) = default;
};

/*
  Define our our population type via sugar template
  In 0.3.1, I introduced the ability to use custom diploid types with the
  library's sugar layer.
*/
using poptype = KTfwd::singlepop<mtype, diploid_t>;

/*
  We will use a gsl_rng_mt19937 as our RNG.
  This type is implicitly convertible to const gsl_rng *,
  and auto-handles the gsl_rng_free steps, etc.
*/
using GSLrng = KTfwd::GSLrng_t<KTfwd::GSL_RNG_MT19937>;

std::size_t
sex_specific_mut_model(std::queue<std::size_t> &mut_recycling_bin,
                       poptype::mcont_t &mutations, const gsl_rng *r,
                       poptype::lookup_table_t &lookup, const double &mu_total,
                       const double &mu_male, const double &mu_female,
                       const double &sigma)
{
    double pos = gsl_rng_uniform(r);
    while (lookup.find(pos) != lookup.end())
        {
            pos = gsl_rng_uniform(r);
        }
    lookup.insert(pos);
    double u = gsl_rng_uniform(r);
    if (u <= mu_male / mu_total)
        {
            return KTfwd::fwdpp_internal::recycle_mutation_helper(
                mut_recycling_bin, mutations, pos, gsl_ran_gaussian(r, sigma),
                false, false);
        }
    else if (u <= (mu_male + mu_female) / mu_total)
        {
            return KTfwd::fwdpp_internal::recycle_mutation_helper(
                mut_recycling_bin, mutations, pos, gsl_ran_gaussian(r, sigma),
                true, false);
        }
    // Otherwise, neutral mutation
    // We "hack" this and assign the mutation a "male" type,
    // As they'll never be used in a fitness calc,
    // as they'll be stored in mutations rather than
    // smutations
    return KTfwd::fwdpp_internal::recycle_mutation_helper(
        mut_recycling_bin, mutations, pos, gsl_ran_gaussian(r, sigma), false,
        true);
}

// Now, we need our own "rules"
struct sexSpecificRules
{
    mutable double wbar, mwbar, fwbar;
    mutable std::vector<double> male_fitnesses, female_fitnesses;
    mutable std::vector<size_t> male_indexes, female_indexes;
    mutable KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr male_lookup,
        female_lookup;
    //! \brief Constructor
    sexSpecificRules()
        : wbar(0.), male_fitnesses(std::vector<double>()),
          female_fitnesses(std::vector<double>()),
          male_indexes(std::vector<size_t>()),
          female_indexes(std::vector<size_t>()),
          male_lookup(KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(nullptr)),
          female_lookup(KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(nullptr))
    {
    }
    //! \brief The "fitness manager"
    template <typename dipcont_t, typename gcont_t, typename mcont_t,
              typename fitness_func>
    void
    w(const dipcont_t &diploids, gcont_t &gametes, mcont_t &mutations,
      const fitness_func &ff) const
    {
        unsigned N_curr = diploids.size();
        if (male_fitnesses.size() < N_curr)
            {
                male_fitnesses.resize(N_curr);
                male_indexes.resize(N_curr);
            }
        if (female_fitnesses.size() < N_curr)
            {
                female_fitnesses.resize(N_curr);
                female_indexes.resize(N_curr);
            }
        wbar = mwbar = fwbar = 0.;

        double w;
        unsigned male = 0, female = 0;
        for (unsigned i = 0; i < N_curr; ++i)
            {
                gametes[diploids[i].first].n = gametes[diploids[i].second].n
                    = 0;
                w = ff(diploids[i], gametes, mutations);
                if (!diploids[i].sex)
                    {
                        male_fitnesses[male] = w;
                        mwbar += w;
                        male_indexes[male++] = i;
                    }
                else
                    {
                        female_fitnesses[female] = w;
                        fwbar += w;
                        female_indexes[female++] = i;
                    }
                wbar += w;
            }
        wbar /= double(diploids.size());
        mwbar /= double(male);
        fwbar /= double(female);

        /*
          Biological point:

          We want to ensure that ANY individual is chosen as a parent
          proportional to w/wbar (conditional on the
          individual being of the desired "sex".

          Thus, we need to transform all male and female fitnesses by wbar.

          If we did not do this transform, individuals would be chosen
          according to w/wbar_sex, which is desired in
          some cases, but not here...
         */
        std::transform(
            male_fitnesses.begin(), male_fitnesses.begin() + male,
            male_fitnesses.begin(),
            std::bind(std::divides<double>(), std::placeholders::_1, wbar));
        std::transform(
            female_fitnesses.begin(), female_fitnesses.begin() + female,
            female_fitnesses.begin(),
            std::bind(std::divides<double>(), std::placeholders::_1, wbar));
        /*!
          Black magic alert:
          fwdpp_internal::gsl_ran_discrete_t_ptr contains a std::unique_ptr
          wrapping the GSL pointer.
          This type has its own deleter, which is convenient, because
          operator= for unique_ptrs automagically calls the deleter before
          assignment!
          Details:
          http://www.cplusplus.com/reference/memory/unique_ptr/operator=
        */
        male_lookup = KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(
            gsl_ran_discrete_preproc(male, &male_fitnesses[0]));
        female_lookup = KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(
            gsl_ran_discrete_preproc(female, &female_fitnesses[0]));
    }

    //! \brief Pick parent one
    inline size_t
    pick1(const gsl_rng *r) const
    {
        return male_indexes[gsl_ran_discrete(r, male_lookup.get())];
    }

    //! \brief Pick parent 2.  Parent 1's data are passed along for models
    //! where that is relevant
    template <typename diploid_t, typename gcont_t, typename mcont_t>
    inline size_t
    pick2(const gsl_rng *r, const size_t &, const double &, const diploid_t &,
          const gcont_t &, const mcont_t &) const
    {
        return female_indexes[gsl_ran_discrete(r, female_lookup.get())];
    }

    //! \brief Update some property of the offspring based on properties of the
    //! parents
    template <typename diploid_t, typename gcont_t, typename mcont_t>
    void
    update(const gsl_rng *r, diploid_t &offspring, const diploid_t &,
           const diploid_t &, const gcont_t &, const mcont_t &) const
    {
        offspring.sex = (gsl_rng_uniform(r) < 0.5);
        return;
    }
};

// We need a fitness model
double
sex_specific_fitness(const poptype::diploid_t &dip,
                     const poptype::gcont_t &gametes,
                     const poptype::mcont_t &mutations, const gsl_rng *r,
                     const double &sigmaE)
{
    double trait_value = std::accumulate(
        gametes[dip.first].smutations.begin(),
        gametes[dip.first].smutations.end(), 0.,
        [&dip, &mutations](const double a, const std::size_t &m) {
            return a + ((dip.sex == mutations[m].sex) ? mutations[m].s : 0.);
        });
    trait_value += std::accumulate(
        gametes[dip.second].smutations.begin(),
        gametes[dip.second].smutations.end(), 0.,
        [&dip, &mutations](const double a, const std::size_t &m) {
            return a + ((dip.sex == mutations[m].sex) ? mutations[m].s : 0.);
        });
    return std::exp(-std::pow(trait_value + gsl_ran_gaussian(r, sigmaE), 2.)
                    / 2.);
}

int
main(int argc, char **argv)
{
    if (argc != 12)
        {
            std::cerr << "Too few arguments.\n"
                      << "Usage: " << argv[0]
                      << " N mu_neutral mu_male mu_female sigma_mu sigma_e "
                         "recrate ngens samplesize nreps seed\n";
            exit(10);
        }
    int argument = 1;

    const unsigned N = atoi(argv[argument++]);
    const double mu_neutral = atof(argv[argument++]);
    const double mu_male = atof(argv[argument++]);
    const double mu_female = atof(argv[argument++]);
    const double sigma = atof(argv[argument++]);
    const double sigmaE = atof(argv[argument++]);
    const double recrate = atof(argv[argument++]);
    const unsigned ngens = atoi(argv[argument++]);
    const unsigned samplesize1 = atoi(argv[argument++]);
    const unsigned nreps = atoi(argv[argument++]);
    const unsigned seed = atoi(argv[argument++]);

    double theta = 4. * (double(N)) * (mu_male + mu_female);
    GSLrng rng(seed);
    std::function<double(void)> recmap
        = std::bind(gsl_rng_uniform, rng.get()); // uniform crossover map
    const double mu_total = mu_neutral + mu_male + mu_female;
    sexSpecificRules rules;
    for (unsigned rep = 0; rep < nreps; ++rep)
        {
            poptype pop(N);
            pop.mutations.reserve(
                size_t(std::ceil(std::log(2 * N) * theta + 0.667 * theta)));
            // Assign "sex"
            for (auto dip = pop.diploids.begin(); dip != pop.diploids.end();
                 ++dip)
                {
                    dip->sex = (gsl_rng_uniform(rng.get())
                                < 0.5); // false = male, true = female.
                }
            for (unsigned generation = 0; generation < ngens; ++generation)
                {
                    // double wbar = KTfwd::sample_diploid(rng,
                    double wbar = KTfwd::experimental::sample_diploid(
                        rng.get(), pop.gametes, pop.diploids, pop.mutations,
                        pop.mcounts, N, mu_total,
                        std::bind(sex_specific_mut_model,
                                  std::placeholders::_1, std::placeholders::_2,
                                  rng.get(), std::ref(pop.mut_lookup),
                                  mu_total, mu_male, mu_female, sigma),
                        std::bind(KTfwd::poisson_xover(), rng.get(), recrate,
                                  0., 1., std::placeholders::_1,
                                  std::placeholders::_2,
                                  std::placeholders::_3),
                        std::bind(sex_specific_fitness, std::placeholders::_1,
                                  std::placeholders::_2, std::placeholders::_3,
                                  rng.get(), sigmaE),
                        pop.neutral, pop.selected,
                        0., // Gotta pass the "selfing" rate, even though it
                        // makes no sense for this model.  API tradeoff for
                        // flexibility...
                        rules);
                    KTfwd::update_mutations(
                        pop.mutations, pop.fixations, pop.fixation_times,
                        pop.mut_lookup, pop.mcounts, generation, 2 * pop.N);
                }

            // Take a sample of size samplesize1.  Two data blocks are
            // returned, one for neutral mutations, and one for selected
            std::pair<std::vector<std::pair<double, std::string>>,
                      std::vector<std::pair<double, std::string>>>
                sample = KTfwd::ms_sample_separate(rng.get(), pop.mutations,
                                                   pop.gametes, pop.diploids,
                                                   samplesize1);

#ifdef HAVE_LIBSEQUENCE
            Sequence::SimData neutral_muts, selected_muts;
            neutral_muts.assign(sample.first.begin(), sample.first.end());
            selected_muts.assign(sample.second.begin(), sample.second.end());
            std::cout << neutral_muts << '\n' << selected_muts << '\n';
#endif
        }
}
