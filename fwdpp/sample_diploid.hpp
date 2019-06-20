#ifndef FWDPP_SAMPLE_DIPLOID_HPP
#define FWDPP_SAMPLE_DIPLOID_HPP

#include <utility>
#include <vector>
#include <fwdpp/fwd_functional.hpp>
#include <fwdpp/insertion_policies.hpp>
namespace fwdpp
{
    /*! \brief Sample the next generation of dipliods in an individual-based
      simulation.  Constant population size case.
      \param r GSL random number generator
      \param haploid_genomes Gametes currently in population
      \param diploids Vector of parents from which we sample offspring
      \param mutations Mutations currently in population
      \param mcounts Vector of integers corresponding to counts of each element
      in mutations
      \param N_curr The population size
      \param mu The total mutation rate per haploid_genome
      \param mmodel Mutation model policy
      \param rec_pol Recombination model policy
      \param ff Policy calculating the fitness of a diploid
      \param neutral
      \param selected
      \param f Probability that a mating is a selfing event
      \param mp Policy determining how whether or not to remove fixed variants
      from the haploid_genomes.

      \note diploids will be updated to reflect the new diploid genotypes
      post-sampling (the descedants).  Gametes will be changed by mutation,
      recombination, and sampling.  Mutations will be changed by mutation and
      sampling.
      \return The mean fitness of the parental generation
      \example diploid_ind.cc
      \example diploid_fixed_sh_ind.cc
    */
    template <typename haploid_genome_type, typename haploid_genome_cont_type_allocator,
              typename mutation_type, typename mutation_cont_type_allocator,
              typename diploid_geno_t, typename diploid_vector_type_allocator,
              typename diploid_fitness_function, typename mutation_model,
              typename recombination_policy,
              template <typename, typename> class haploid_genome_cont_type,
              template <typename, typename> class mutation_cont_type,
              template <typename, typename> class diploid_vector_type,
              typename mutation_removal_policy = std::true_type>
    double sample_diploid(
        const gsl_rng *r,
        haploid_genome_cont_type<haploid_genome_type, haploid_genome_cont_type_allocator> &haploid_genomes,
        diploid_vector_type<diploid_geno_t, diploid_vector_type_allocator>
            &diploids,
        mutation_cont_type<mutation_type, mutation_cont_type_allocator>
            &mutations,
        std::vector<uint_t> &mcounts, const uint_t &N_curr, const double &mu,
        const mutation_model &mmodel, const recombination_policy &rec_pol,
        const diploid_fitness_function &ff,
        typename haploid_genome_type::mutation_container &neutral,
        typename haploid_genome_type::mutation_container &selected,
        const double f = 0.,
        const mutation_removal_policy mp = mutation_removal_policy());

    /*! \brief Sample the next generation of dipliods in an individual-based
      simulation.  Changing population size case.
      \param r GSL random number generator
      \param haploid_genomes Gametes currently in population
      \param diploids Vector of parents from which we sample offspring
      \param mutations Mutations currently in population
      \param mcounts Vector of integers corresponding to counts of each element
      in mutations
      \param N_curr The current population size
      \param N_next The population size after sampling
      \param mu The total mutation rate per haploid_genome
      \param mmodel Mutation model policy
      \param rec_pol Recombination model policy
      \param ff Policy calculating the fitness of a diploid
      \param neutral
      \param selected
      \param f Probability that a mating is a selfing event
      \param mp Policy determining how whether or not to remove fixed variants
      from the haploid_genomes.

      \note diploids will be updated to reflect the new diploid genotypes
      post-sampling (the descedants).  Gametes will be changed by mutation,
      recombination, and sampling.  Mutations will be changed by mutation and
      sampling.
      \return The mean fitness of the parental generation
    */
    template <typename haploid_genome_type, typename haploid_genome_cont_type_allocator,
              typename mutation_type, typename mutation_cont_type_allocator,
              typename diploid_geno_t, typename diploid_vector_type_allocator,
              typename diploid_fitness_function, typename mutation_model,
              typename recombination_policy,
              template <typename, typename> class haploid_genome_cont_type,
              template <typename, typename> class mutation_cont_type,
              template <typename, typename> class diploid_vector_type,
              typename mutation_removal_policy = std::true_type>
    double sample_diploid(
        const gsl_rng *r,
        haploid_genome_cont_type<haploid_genome_type, haploid_genome_cont_type_allocator> &haploid_genomes,
        diploid_vector_type<diploid_geno_t, diploid_vector_type_allocator>
            &diploids,
        mutation_cont_type<mutation_type, mutation_cont_type_allocator>
            &mutations,
        std::vector<uint_t> &mcounts, const uint_t &N_curr,
        const uint_t &N_next, const double &mu, const mutation_model &mmodel,
        const recombination_policy &rec_pol,
        const diploid_fitness_function &ff,
        typename haploid_genome_type::mutation_container &neutral,
        typename haploid_genome_type::mutation_container &selected,
        const double f = 0.,
        const mutation_removal_policy mp = mutation_removal_policy());
} // namespace fwdpp

#include <fwdpp/sample_diploid.tcc>
#endif
