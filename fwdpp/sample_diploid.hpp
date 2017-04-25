#ifndef __FWDPP_SAMPLE_DIPLOID_HPP__
#define __FWDPP_SAMPLE_DIPLOID_HPP__

#include <utility>
#include <vector>
#include <fwdpp/fwd_functional.hpp>
#include <fwdpp/insertion_policies.hpp>
namespace KTfwd
{
    /*! \brief Sample the next generation of dipliods in an individual-based
      simulation.  Constant population size case.
      \param r GSL random number generator
      \param gametes Gametes currently in population
      \param diploids Vector of parents from which we sample offspring
      \param mutations Mutations currently in population
      \param mcounts Vector of integers corresponding to counts of each element
      in mutations
      \param N_curr The population size
      \param mu The total mutation rate per gamete
      \param mmodel Mutation model policy
      \param rec_pol Recombination model policy
      \param ff Policy calculating the fitness of a diploid
      \param neutral
      \param selected
      \param f Probability that a mating is a selfing event
      \param mp Policy determining how whether or not to remove fixed variants
      from the gametes.
      \param gpolicy_mut Policy determining how new gametes are added to
      population after a mutation event

      \note diploids will be updated to reflect the new diploid genotypes
      post-sampling (the descedants).  Gametes will be changed by mutation,
      recombination, and sampling.  Mutations will be changed by mutation and
      sampling.
      \return The mean fitness of the parental generation
      \example diploid_ind.cc
      \example pfix.cc
      \example diploid_fixed_sh_ind_lambda.cc
    */
    template <typename gamete_type, typename gamete_cont_type_allocator,
              typename mutation_type, typename mutation_cont_type_allocator,
              typename diploid_geno_t, typename diploid_vector_type_allocator,
              typename diploid_fitness_function, typename mutation_model,
              typename recombination_policy,
              template <typename, typename> class gamete_cont_type,
              template <typename, typename> class mutation_cont_type,
              template <typename, typename> class diploid_vector_type,
              typename mutation_removal_policy = std::true_type,
              typename gamete_insertion_policy = emplace_back>
    double sample_diploid(
        const gsl_rng *r,
        gamete_cont_type<gamete_type, gamete_cont_type_allocator> &gametes,
        diploid_vector_type<diploid_geno_t, diploid_vector_type_allocator>
            &diploids,
        mutation_cont_type<mutation_type, mutation_cont_type_allocator>
            &mutations,
        std::vector<uint_t> &mcounts, const uint_t &N_curr, const double &mu,
        const mutation_model &mmodel, const recombination_policy &rec_pol,
        const diploid_fitness_function &ff,
        typename gamete_type::mutation_container &neutral,
        typename gamete_type::mutation_container &selected,
        const double f = 0.,
        const mutation_removal_policy mp = mutation_removal_policy(),
        const gamete_insertion_policy &gpolicy_mut
        = gamete_insertion_policy());

    /*! \brief Sample the next generation of dipliods in an individual-based
      simulation.  Changing population size case.
      \param r GSL random number generator
      \param gametes Gametes currently in population
      \param diploids Vector of parents from which we sample offspring
      \param mutations Mutations currently in population
      \param mcounts Vector of integers corresponding to counts of each element
      in mutations
      \param N_curr The current population size
      \param N_next The population size after sampling
      \param mu The total mutation rate per gamete
      \param mmodel Mutation model policy
      \param rec_pol Recombination model policy
      \param ff Policy calculating the fitness of a diploid
      \param neutral
      \param selected
      \param f Probability that a mating is a selfing event
      \param mp Policy determining how whether or not to remove fixed variants
      from the gametes.
      \param gpolicy_mut Policy determining how new gametes are added to
      population after a mutation event

      \note diploids will be updated to reflect the new diploid genotypes
      post-sampling (the descedants).  Gametes will be changed by mutation,
      recombination, and sampling.  Mutations will be changed by mutation and
      sampling.
      \return The mean fitness of the parental generation
      \example bneck_selection_ind.cc
    */
    template <typename gamete_type, typename gamete_cont_type_allocator,
              typename mutation_type, typename mutation_cont_type_allocator,
              typename diploid_geno_t, typename diploid_vector_type_allocator,
              typename diploid_fitness_function, typename mutation_model,
              typename recombination_policy,
              template <typename, typename> class gamete_cont_type,
              template <typename, typename> class mutation_cont_type,
              template <typename, typename> class diploid_vector_type,
              typename mutation_removal_policy = std::true_type,
              typename gamete_insertion_policy = emplace_back>
    double sample_diploid(
        const gsl_rng *r,
        gamete_cont_type<gamete_type, gamete_cont_type_allocator> &gametes,
        diploid_vector_type<diploid_geno_t, diploid_vector_type_allocator>
            &diploids,
        mutation_cont_type<mutation_type, mutation_cont_type_allocator>
            &mutations,
        std::vector<uint_t> &mcounts, const uint_t &N_curr,
        const uint_t &N_next, const double &mu, const mutation_model &mmodel,
        const recombination_policy &rec_pol,
        const diploid_fitness_function &ff,
        typename gamete_type::mutation_container &neutral,
        typename gamete_type::mutation_container &selected,
        const double f = 0.,
        const mutation_removal_policy mp = mutation_removal_policy(),
        const gamete_insertion_policy &gpolicy_mut
        = gamete_insertion_policy());

    /*! \brief Evolve a metapopulation where demes are not changing size.  For
      individual-based sims.
      Evolve a metapopulation where demes are not changing size.  For
      individual-based sims.
      \param r GSL random number generator
      \param gametes Container of gametes
      \param diploids Vector of vectors of parents from which we sample
      offspring
      \param mutations Container mutations currently in population
      \param mcounts Vector of integers corresponding to counts of each element
      in mutations
      \param N_curr The population sizes.  There must be diploids.size() values
      in this array
      \param mu The total mutation rate per gamete
      \param mmodel Mutation model policy
      \param rec_pol Recombination model policy
      \param gpolicy_mut Policy determining how new gametes are added to
      population after a mutation event
      \param ffs Container of fitness functions.  One for each deme.
      \param neutral
      \param selected
      \param mp Policy determining how to remove mutations from a diploid
      (e.g., removing fixed and/or lost mutations)
      \param mig Migration policy.  This function/function object must take a
      single size_t (values 0 to diploids.size()-1).  If no migration event
      occurs, the passed value is returned.  Otherwise, a size_t representing
      the index of the deme from which the other parent comes (aka the migrant)
      is returned.
      \param f Probability that a mating is a selfing event.  This is an array,
      with 1 f per deme.

      \note diploids will be updated to reflect the new diploid genotypes
      post-sampling (the descedants).  Gametes will be changed by mutation,
      recombination, and sampling.  Mutations will be changed by mutation and
      sampling.
      \return The mean fitness of the parental generation

      \example migsel_ind.cc
    */
    template <typename gamete_type, typename mutation_type,
              typename metapop_diploid_vector_type_allocator,
              typename gamete_cont_type_allocator,
              typename mutation_cont_type_allocator, typename diploid_geno_t,
              typename diploid_vector_type_allocator,
              typename diploid_fitness_function_container,
              typename mutation_model, typename recombination_policy,
              typename migration_policy,
              template <typename, typename> class gamete_cont_type,
              template <typename, typename> class mutation_cont_type,
              template <typename, typename> class diploid_vector_type,
              template <typename, typename> class metapop_diploid_vector_type,
              typename mutation_removal_policy = std::true_type,
              typename gamete_insertion_policy = emplace_back>
    std::vector<double> sample_diploid(
        const gsl_rng *r,
        gamete_cont_type<gamete_type, gamete_cont_type_allocator> &gametes,
        metapop_diploid_vector_type<diploid_vector_type<diploid_geno_t,
                                                        diploid_vector_type_allocator>,
                                    metapop_diploid_vector_type_allocator>
            &diploids,
        mutation_cont_type<mutation_type, mutation_cont_type_allocator>
            &mutations,
        std::vector<uint_t> &mcounts, const uint_t *N_curr, const double &mu,
        const mutation_model &mmodel, const recombination_policy &rec_pol,
        const diploid_fitness_function_container &ffs,
        const migration_policy &mig,
        typename gamete_type::mutation_container &neutral,
        typename gamete_type::mutation_container &selected,
        const double *f = nullptr,
        const mutation_removal_policy &mp = mutation_removal_policy(),
        const gamete_insertion_policy &gpolicy_mut
        = gamete_insertion_policy());

    /*! \brief Evolve a metapopulation where demes may be changing size.  For
      individual-based sims.
      Evolve a metapopulation where demes may be changing size.  For
      individual-based sims.
      \param r GSL random number generator
      \param gametes Container of gametes
      \param diploids Vector of vectors of parents from which we sample
      offspring
      \param mutations Container mutations currently in population
      \param mcounts Vector of integers corresponding to counts of each element
      in mutations
      \param N_curr The current population sizes.  There must be
      diploids.size() values in this array
      \param N_next The population sizes in the daughter generation.  There
      must be diploids->size() values in this array
      \param mu The total mutation rate per gamete
      \param mmodel Mutation model policy
      \param rec_pol Recombination model policy
      \param gpolicy_mut Policy determining how new gametes are added to
      population after a mutation event
      \param ffs Container of fitness functions.  One for each deme.
      \param neutral
      \param selected
      \param mp Policy determining how to remove mutations from a diploid
      (e.g., removing fixed and/or lost mutations)
      \param mig Migration policy.  This function/function object must take a
      single size_t (values 0 to metapop->size()-1).  If no migration event
      occurs, the passed value is returned.  Otherwise, a size_t representing
      the index of the deme from which the other parent comes (aka the migrant)
      is returned.
      \param f Probability that a mating is a selfing event.  This is an array,
      with 1 f per deme.

      \note diploids will be updated to reflect the new diploid genotypes
      post-sampling (the descedants).  Gametes will be changed by mutation,
      recombination, and sampling.  Mutations will be changed by mutation and
      sampling.
      \return The mean fitness of the parental generation
    */
    template <typename gamete_type, typename mutation_type,
              typename metapop_diploid_vector_type_allocator,
              typename gamete_cont_type_allocator,
              typename mutation_cont_type_allocator, typename diploid_geno_t,
              typename diploid_vector_type_allocator,
              typename diploid_fitness_function_container,
              typename mutation_model, typename recombination_policy,
              typename migration_policy,
              template <typename, typename> class gamete_cont_type,
              template <typename, typename> class mutation_cont_type,
              template <typename, typename> class diploid_vector_type,
              template <typename, typename> class metapop_diploid_vector_type,
              typename mutation_removal_policy = std::true_type,
              typename gamete_insertion_policy = emplace_back>
    std::vector<double> sample_diploid(
        const gsl_rng *r,
        gamete_cont_type<gamete_type, gamete_cont_type_allocator> &gametes,
        metapop_diploid_vector_type<diploid_vector_type<diploid_geno_t,
                                                        diploid_vector_type_allocator>,
                                    metapop_diploid_vector_type_allocator>
            &diploids,
        mutation_cont_type<mutation_type, mutation_cont_type_allocator>
            &mutations,
        std::vector<uint_t> &mcounts, const uint_t *N_curr,
        const uint_t *N_next, const double &mu, const mutation_model &mmodel,
        const recombination_policy &rec_pol,
        const diploid_fitness_function_container &ffs,
        const migration_policy &mig,
        typename gamete_type::mutation_container &neutral,
        typename gamete_type::mutation_container &selected,
        const double *f = nullptr,
        const mutation_removal_policy &mp = mutation_removal_policy(),
        const gamete_insertion_policy &gpolicy_mut
        = gamete_insertion_policy());

    /*! \brief Single deme, multilocus model, changing population size
     */
    template <
        typename diploid_geno_t, typename gamete_type,
        typename gamete_cont_type_allocator, typename mutation_type,
        typename mutation_cont_type_allocator,
        typename diploid_vector_type_allocator,
        typename locus_vector_type_allocator,
        typename diploid_fitness_function, typename mutation_model_container,
        typename recombination_policy_container,
        template <typename, typename> class gamete_cont_type,
        template <typename, typename> class mutation_cont_type,
        template <typename, typename> class diploid_vector_type,
        template <typename, typename> class locus_vector_type,
        typename mutation_removal_policy = std::true_type,
        typename gamete_insertion_policy = emplace_back>
    double sample_diploid(
        const gsl_rng *r,
        gamete_cont_type<gamete_type, gamete_cont_type_allocator> &gametes,
        diploid_vector_type<locus_vector_type<diploid_geno_t,
                                              locus_vector_type_allocator>,
                            diploid_vector_type_allocator> &diploids,
        mutation_cont_type<mutation_type, mutation_cont_type_allocator>
            &mutations,
        std::vector<uint_t> &mcounts, const uint_t &N_curr,
        const uint_t &N_next, const double *mu,
        const mutation_model_container &mmodel,
        const recombination_policy_container &rec_policies,
		const std::vector<std::function<unsigned(void)>> & interlocus_rec,
        const diploid_fitness_function &ff,
        typename gamete_type::mutation_container &neutral,
        typename gamete_type::mutation_container &selected,
        const double &f = 0,
        const mutation_removal_policy &mp = mutation_removal_policy(),
        const gamete_insertion_policy &gpolicy_mut
        = gamete_insertion_policy());

    /*! \brief Single deme, multilocus model, constant population size
      \example diploid_ind_2locus.cc
    */
    // single deme, constant N
    template <
        typename diploid_geno_t, typename gamete_type,
        typename gamete_cont_type_allocator, typename mutation_type,
        typename mutation_cont_type_allocator,
        typename diploid_vector_type_allocator,
        typename locus_vector_type_allocator,
        typename diploid_fitness_function, typename mutation_model_container,
        typename recombination_policy_container,
        template <typename, typename> class gamete_cont_type,
        template <typename, typename> class mutation_cont_type,
        template <typename, typename> class diploid_vector_type,
        template <typename, typename> class locus_vector_type,
        typename mutation_removal_policy = std::true_type,
        typename gamete_insertion_policy = emplace_back>
    double sample_diploid(
        const gsl_rng *r,
        gamete_cont_type<gamete_type, gamete_cont_type_allocator> &gametes,
        diploid_vector_type<locus_vector_type<diploid_geno_t,
                                              locus_vector_type_allocator>,
                            diploid_vector_type_allocator> &diploids,
        mutation_cont_type<mutation_type, mutation_cont_type_allocator>
            &mutations,
        std::vector<uint_t> &mcounts, const uint_t &N, const double *mu,
        const mutation_model_container &mmodel,
        const recombination_policy_container &rec_policies,
		const std::vector<std::function<unsigned(void)>> & interlocus_rec,
        const diploid_fitness_function &ff,
        typename gamete_type::mutation_container &neutral,
        typename gamete_type::mutation_container &selected,
        const double &f = 0,
        const mutation_removal_policy &mp = mutation_removal_policy(),
        const gamete_insertion_policy &gpolicy_mut
        = gamete_insertion_policy());
}

#include <fwdpp/sample_diploid.tcc>
#endif
