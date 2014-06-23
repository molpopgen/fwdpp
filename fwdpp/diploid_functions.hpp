#ifndef _DIPLOID_FUNCTIONS_HPP_
#define _DIPLOID_FUNCTIONS_HPP_

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/*! \file diploid_functions.hpp
  \brief Wright-Fisher sampling in a finite population and recombination 
*/ 

namespace KTfwd
{
  /*! \brief Wright-Fisher sampling of a single population of constant size
    Performs multinomial sampling of gametes proportional to marginal fitness divided by 
    population mean fitness under constant population size

    \param r GSL random number generator
    \param gametes Vector of gametes in the population
    \param twoN Twice the diploid population size
    \param ff A fitness function taking two gamete_type as arguments
    \param mp A policy determining how mutations are removed from populatio after sampling.  See KTfwd::mutation_remover for example
    \param f  Probability that a mating is an inbreeding event

    \example diploid.cc
    \example RHH.cc
  */
  template< typename gamete_type,
	    typename vector_type_allocator,
	    template<typename,typename> class vector_type,
	    typename diploid_fitness_function,
	    typename mutation_removal_policy>
  double sample_diploid(gsl_rng * r, vector_type<gamete_type,vector_type_allocator > * gametes, 
			const unsigned & twoN, 
			const diploid_fitness_function & ff,
			const mutation_removal_policy & mp,
			const double & f = 0.);
  
  /*! \brief Wright-Fisher sampling of a single population of changing size
    Performs multinomial sampling of gametes proportional to marginal fitness divided by 
    population mean fitness where population size changes from twoN_curr to twoN_next
    
    \param r GSL random number generator
    \param gametes Vector of gametes in the population
    \param twoN_curr Twice the diploid population size in the current generation
    \param twoN_next Twice the diploid population size in the next generation
    \param ff A fitness function taking two gamete_type as arguments
    \param mp A policy determining how mutations are removed from populatio after sampling.  See KTfwd::mutation_remover for example
    \param f  Probability that a mating is an inbreeding event
  */
  template< typename gamete_type,
	    typename vector_type_allocator,
	    template<typename,typename> class vector_type,
	    typename diploid_fitness_function,
	    typename mutation_removal_policy>
  double sample_diploid(gsl_rng * r, vector_type<gamete_type,vector_type_allocator > * gametes, 
			const unsigned & twoN_curr, const unsigned & twoN_next,
			const diploid_fitness_function & ff,
			const mutation_removal_policy & mp,
			const double & f);

/*! \brief Wright-Fisher sampling of a metapopulation population of constant size
  Performs multinomial sampling of gametes proportional to marginal fitness divided by 
  population mean fitness under constant population size for a metapopulation
  
  \param r GSL random number generator
  \param metapop Vector vector of gametes in the metapopulation
  \param twoNs Array of 2N for each deme
  \param metapopsize Sum of values in [twoNs,twoNs+metapop->size())
  \param ffs A vector fitness functions taking two gamete_type as arguments.  Each deme can have a different fitness function
  \param mp A policy determining how mutations are removed from populatio after sampling.  See KTfwd::mutation_remover for example
  \example diploid_twopop_mig.cc
*/
template< typename gamete_type,
	  typename vector_type_allocator1,
	  typename vector_type_allocator2,
	  template<typename,typename> class vector_type1,
	  template<typename,typename> class vector_type2,
	  typename diploid_fitness_function_container,
	  typename mutation_removal_policy>
std::vector<double> sample_diploid(gsl_rng * r, 
				   vector_type2< vector_type1<gamete_type, vector_type_allocator1 >, vector_type_allocator2 > * metapop,
				   const unsigned * twoNs,
				   const unsigned & metapopsize,
				   const diploid_fitness_function_container & ffs,
				   const mutation_removal_policy & mp);
  
  /*!  \brief Wright-Fisher sampling of a metapopulation population of constant size with inbreeding
    Performs multinomial sampling of gametes proportional to marginal fitness divided by 
    population mean fitness under constant population size for a metapopulation
    
    \param r GSL random number generator
    \param metapop Vector vector of gametes in the metapopulation
    \param twoNs Array of 2N for each deme in the current generation
    \param twoNs_next Array of 2N for each deme in the next generation
    \param metapopsize Sum of values in [twoNs,twoNs+metapop->size())
    \param ffs A vector fitness functions taking two gamete_type as arguments.  Each deme can have a different fitness function
    \param mp A policy determining how mutations are removed from populatio after sampling.  See KTfwd::mutation_remover for example
  */
  template< typename gamete_type,
	    typename vector_type_allocator1,
	    typename vector_type_allocator2,
	    template<typename,typename> class vector_type1,
	    template<typename,typename> class vector_type2,
	    typename diploid_fitness_function_container,
	    typename mutation_removal_policy>
  std::vector<double> sample_diploid(gsl_rng * r, 
				     vector_type2< vector_type1<gamete_type, vector_type_allocator1 >, vector_type_allocator2 > * metapop,
				     const unsigned * twoNs,
				     const unsigned * twoNs_next,
				     const unsigned & metapopsize,
				     const diploid_fitness_function_container & ffs,
				     const mutation_removal_policy & mp);

/*! \brief Wright-Fisher sampling of a metapopulation population of changing size
  Performs multinomial sampling of gametes proportional to marginal fitness divided by 
  population mean fitness under constant population size for a metapopulation
  
  \param r GSL random number generator
  \param metapop Vector vector of gametes in the metapopulation
  \param twoNs Array of 2N for each deme in the current generation
  \param metapopsize Sum of values in [twoNs,twoNs+metapop->size())
  \param ffs A vector fitness functions taking two gamete_type as arguments.  Each deme can have a different fitness function
  \param mp A policy determining how mutations are removed from populatio after sampling.  See KTfwd::mutation_remover for example
  \param fs  An array of inbreeding probabilities for each deme.
*/
template< typename gamete_type,
	  typename vector_type_allocator1,
	  typename vector_type_allocator2,
	  template<typename,typename> class vector_type1,
	  template<typename,typename> class vector_type2,
	  typename diploid_fitness_function_container,
	  typename mutation_removal_policy>
std::vector<double> sample_diploid(gsl_rng * r, 
				   vector_type2< vector_type1<gamete_type, vector_type_allocator1 >, vector_type_allocator2 > * metapop,
				   const unsigned * twoNs,
				   const unsigned & metapopsize,
				   const diploid_fitness_function_container & ffs,
				   const mutation_removal_policy & mp,
				   const double * fs);
  
  /*! \brief Wright-Fisher sampling of a metapopulation population of changing size with inbreeding
    Performs multinomial sampling of gametes proportional to marginal fitness divided by 
    population mean fitness under constant population size for a metapopulation
    
    \param r GSL random number generator
    \param metapop Vector vector of gametes in the metapopulation
    \param twoNs Array of 2N for each deme in the current generation
    \param twoNs_next Array of 2N for each deme in the next generation
    \param metapopsize Sum of values in [twoNs,twoNs+metapop->size())
    \param ffs A vector fitness functions taking two gamete_type as arguments.  Each deme can have a different fitness function
    \param mp A policy determining how mutations are removed from populatio after sampling.  See KTfwd::mutation_remover for example
    \param fs  An array of inbreeding probabilities for each deme.
  */
  template< typename gamete_type,
	    typename vector_type_allocator1,
	    typename vector_type_allocator2,
	    template<typename,typename> class vector_type1,
	    template<typename,typename> class vector_type2,
	    typename diploid_fitness_function_container,
	    typename mutation_removal_policy>
  std::vector<double> sample_diploid(gsl_rng * r, 
				     vector_type2< vector_type1<gamete_type, vector_type_allocator1 >, vector_type_allocator2 > * metapop,
				     const unsigned * twoNs,
				     const unsigned * twoNs_next,
				     const unsigned & metapopsize,
				     const diploid_fitness_function_container & ffs,
				     const mutation_removal_policy & mp,
				     const double * fs);

/*! \brief Recombination
  Recombination
  
  \param r GSL random number generator
  \param gametes Vector of gametes
  \param twoN Twice the number of diploids
  \param littler  The recombination rate PER GAMETE!!! (I.e., 1/2 the usual definition of the recombination rate)
  \param mf A function/function object returning a position along the genetic map
  \param f  Probability of inbreeding
  
  \example diploid.cc
*/
template< typename gamete_type,
	  typename vector_type_allocator,
	  template<typename,typename> class vector_type,
	  typename map_function>
unsigned recombine(gsl_rng * r,
		   vector_type<gamete_type,vector_type_allocator > * gametes,
		   const unsigned & twoN, 
		   const double & littler,
		   const map_function & mf,
		   const double & f = 0);
  
  //PROTOTYPES FOR INDIVIDUAL-BASED SIMULATION ROUTINES ARE BELOW
  
  /*! \brief Sample the next generation of dipliods in an individual-based simulation.  Constant population size case.
     \param r GSL random number generator
     \param gametes Pointer to list of gametes currently in population
     \param diploids Pointer to the vector of parents from which we sample offspring
     \param mutations Pointer to list of mutations currently in population
     \param N_curr The population size
     \param mu The total mutation rate per gamete
     \param mmodel Mutation model policy
     \param rec_pol Recombination model policy
     \param mpolicy Policy determining how new mutations are added to the population
     \param gpolicy_mut Policy determining how new gametes are added to population after a mutation event
     \param ff Policy calculating the fitness of a diploid
     \param mp Policy determining how to remove mutations from a diploid (e.g., removing fixed and/or lost mutations)
     \param f Probability that a mating is a selfing event

     \note diploids will be updated to reflect the new diploid genotypes post-sampling (the descedants).  Gametes will be changed by mutation, recombination, and sampling.  Mutations will be changed by mutation and sampling.
     \return The mean fitness of the parental generation
     \example diploid_ind.cc
     \example pfix.cc
   */
template< typename gamete_type,
	  typename gamete_list_type_allocator,
	  typename mutation_list_type_allocator,
	  typename diploid_vector_type_allocator,
	  typename diploid_fitness_function,
	  typename mutation_removal_policy,
	  typename mutation_model,
	  typename recombination_policy,
	  typename mutation_insertion_policy,
	  typename gamete_insertion_policy,
	  template<typename,typename> class gamete_list_type,
	  template<typename,typename> class mutation_list_type,
	  template<typename,typename> class diploid_vector_type>
  double
  sample_diploid(gsl_rng * r,
		 gamete_list_type<gamete_type,gamete_list_type_allocator > * gametes,
		 diploid_vector_type<std::pair<typename gamete_list_type< gamete_type,gamete_list_type_allocator >::iterator,
					       typename gamete_list_type< gamete_type,gamete_list_type_allocator >::iterator>,
		 diploid_vector_type_allocator> * diploids,
		 mutation_list_type<typename gamete_type::mutation_type,mutation_list_type_allocator > * mutations, 
		 const unsigned & N_curr, 
		 const double & mu,
		 const mutation_model & mmodel,
		 const recombination_policy & rec_pol,
		 const mutation_insertion_policy & mpolicy,
		 const gamete_insertion_policy & gpolicy_mut,
		 const diploid_fitness_function & ff,
		 const mutation_removal_policy & mp,
		 const double & f = 0);
 
  /*! \brief Sample the next generation of dipliods in an individual-based simulation.  Changing population size case.
    \param r GSL random number generator
    \param gametes Pointer to list of gametes currently in population
    \param diploids Pointer to the vector of parents from which we sample offspring
    \param mutations Pointer to list of mutations currently in population
    \param N_curr The population size in the current generation
    \param N_next The population size in the next generation
    \param mu The total mutation rate per gamete
    \param mmodel Mutation model policy
    \param rec_pol Recombination model policy
    \param mpolicy Policy determining how new mutations are added to the population
    \param gpolicy_mut Policy determining how new gametes are added to population after a mutation event
    \param ff Policy calculating the fitness of a diploid
    \param mp Policy determining how to remove mutations from a diploid (e.g., removing fixed and/or lost mutations)
    \param f Probability that a mating is a selfing event
    
     \note diploids will be updated to reflect the new diploid genotypes post-sampling (the descedants).  Gametes will be changed by mutation, recombination, and sampling.  Mutations will be changed by mutation and sampling.
     \return The mean fitness of the parental generation
     \example bneck_selection.cc 
     \example TFL2013.cc
   */
template< typename gamete_type,
	  typename gamete_list_type_allocator,
	  typename mutation_list_type_allocator,
	  typename diploid_vector_type_allocator,
	  typename diploid_fitness_function,
	  typename mutation_removal_policy,
	  typename mutation_model,
	  typename recombination_policy,
	  typename mutation_insertion_policy,
	  typename gamete_insertion_policy,
	  template<typename,typename> class gamete_list_type,
	  template<typename,typename> class mutation_list_type,
	  template<typename,typename> class diploid_vector_type>
  double
  sample_diploid(gsl_rng * r,
		 gamete_list_type<gamete_type,gamete_list_type_allocator > * gametes,
		 diploid_vector_type<std::pair<typename gamete_list_type< gamete_type,gamete_list_type_allocator >::iterator,
					       typename gamete_list_type< gamete_type,gamete_list_type_allocator >::iterator>,
		 diploid_vector_type_allocator> * diploids,
		 mutation_list_type<typename gamete_type::mutation_type,mutation_list_type_allocator > * mutations, 
		 const unsigned & N_curr, 
		 const unsigned & N_next, 
		 const double & mu,
		 const mutation_model & mmodel,
		 const recombination_policy & rec_pol,
		 const mutation_insertion_policy & mpolicy,
		 const gamete_insertion_policy & gpolicy_mut,
		 const diploid_fitness_function & ff,
		 const mutation_removal_policy & mp,
		 const double & f = 0);

  /*! \brief Evolve a metapopulation where demes are not changing size.  For individual-based sims.
    Evolve a metapopulation where demes are not changing size.  For individual-based sims.
     \param r GSL random number generator
     \param metapop Pointer to vector of lists of gametes.
     \param dipliods Pointer to the vector of vectors of parents from which we sample offspring
     \param mutations Pointer to list of mutations currently in population
     \param N_curr The population sizes.  There must be diploids->size() values in this array
     \param mu The total mutation rate per gamete
     \param mmodel Mutation model policy
     \param rec_pol Recombination model policy
     \param mpolicy Policy determining how new mutations are added to the population
     \param gpolicy_mut Policy determining how new gametes are added to population after a mutation event
     \param ffs Container of fitness functions.  One for each deme.
     \param mp Policy determining how to remove mutations from a diploid (e.g., removing fixed and/or lost mutations)
     \param mig Migration policy.  This function/function object must take a single size_t (values 0 to metapop->size()-1).  If no migration event occurs, the passed value is returned.  Otherwise, a size_t representing the index of the deme from which the other parent comes (aka the migrant) is returned.
     \param f Probability that a mating is a selfing event.  This is an array, with 1 f per deme.

     \note diploids will be updated to reflect the new diploid genotypes post-sampling (the descedants).  Gametes will be changed by mutation, recombination, and sampling.  Mutations will be changed by mutation and sampling.
     \return The mean fitness of the parental generation

     \example migsel_ind.cc
   */
template< typename gamete_type,
	  typename metapop_gamete_vector_type_allocator,
	  typename metapop_diploid_vector_type_allocator,
	  typename gamete_list_type_allocator,
	  typename mutation_list_type_allocator,
	  typename diploid_vector_type_allocator,
	  typename diploid_fitness_function_container,
	  typename mutation_removal_policy,
	  typename mutation_model,
	  typename recombination_policy,
	  typename migration_policy,
	  typename mutation_insertion_policy,
	  typename gamete_insertion_policy,
	  template<typename,typename> class gamete_list_type,
	  template<typename,typename> class mutation_list_type,
	  template<typename,typename> class diploid_vector_type,
	  template<typename,typename> class metapop_gamete_vector_type,
	  template<typename,typename> class metapop_diploid_vector_type>
std::vector< double >
sample_diploid(gsl_rng * r,
	       metapop_gamete_vector_type < gamete_list_type<gamete_type,gamete_list_type_allocator > ,
	       metapop_gamete_vector_type_allocator > * metapop,
  metapop_diploid_vector_type < diploid_vector_type<std::pair<typename gamete_list_type< gamete_type,gamete_list_type_allocator >::iterator,
							      typename gamete_list_type< gamete_type,gamete_list_type_allocator >::iterator>,
						    diploid_vector_type_allocator>,
				metapop_diploid_vector_type_allocator > * diploids,
	       mutation_list_type<typename gamete_type::mutation_type,mutation_list_type_allocator > * mutations, 
	       const unsigned * N_curr, 
	       const double & mu,
	       const mutation_model & mmodel,
	       const recombination_policy & rec_pol,
	       const mutation_insertion_policy & mpolicy,
	       const gamete_insertion_policy & gpolicy_mut,
	       const diploid_fitness_function_container & ffs,
	       const mutation_removal_policy & mp,
	       const migration_policy & mig,
	       const double * f = NULL);

  /*! \brief Evolve a metapopulation where demes may be changing size.  For individual-based sims.
    Evolve a metapopulation where demes may be changing size.  For individual-based sims.
    \param r GSL random number generator
    \param metapop Pointer to vector of lists of gametes.
    \param diploids Pointer to the vector of vectors of parents from which we sample offspring
    \param mutations Pointer to list of mutations currently in population
    \param N_curr The current population sizes.  There must be diploids->size() values in this array
    \param N_next The population sizes in the daughter generation.  There must be diploids->size() values in this array
    \param mu The total mutation rate per gamete
    \param mmodel Mutation model policy
    \param rec_pol Recombination model policy
    \param mpolicy Policy determining how new mutations are added to the population
    \param gpolicy_mut Policy determining how new gametes are added to population after a mutation event
    \param ffs Container of fitness functions.  One for each deme.
    \param mp Policy determining how to remove mutations from a diploid (e.g., removing fixed and/or lost mutations)
    \param mig Migration policy.  This function/function object must take a single size_t (values 0 to metapop->size()-1).  If no migration event occurs, the passed value is returned.  Otherwise, a size_t representing the index of the deme from which the other parent comes (aka the migrant) is returned.
    \param f Probability that a mating is a selfing event.  This is an array, with 1 f per deme.
    
    \note diploids will be updated to reflect the new diploid genotypes post-sampling (the descedants).  Gametes will be changed by mutation, recombination, and sampling.  Mutations will be changed by mutation and sampling.
    \return The mean fitness of the parental generation
  */
template< typename gamete_type,
	  typename metapop_gamete_vector_type_allocator,
	  typename metapop_diploid_vector_type_allocator,
	  typename gamete_list_type_allocator,
	  typename mutation_list_type_allocator,
	  typename diploid_vector_type_allocator,
	  typename diploid_fitness_function_container,
	  typename mutation_removal_policy,
	  typename mutation_model,
	  typename recombination_policy,
	  typename migration_policy,
	  typename mutation_insertion_policy,
	  typename gamete_insertion_policy,
	  template<typename,typename> class gamete_list_type,
	  template<typename,typename> class mutation_list_type,
	  template<typename,typename> class diploid_vector_type,
	  template<typename,typename> class metapop_gamete_vector_type,
	  template<typename,typename> class metapop_diploid_vector_type>
std::vector< double >
sample_diploid(gsl_rng * r,
	       metapop_gamete_vector_type < gamete_list_type<gamete_type,gamete_list_type_allocator > ,
	       metapop_gamete_vector_type_allocator > * metapop,
  metapop_diploid_vector_type < diploid_vector_type<std::pair<typename gamete_list_type< gamete_type,gamete_list_type_allocator >::iterator,
							      typename gamete_list_type< gamete_type,gamete_list_type_allocator >::iterator>,
						    diploid_vector_type_allocator>,
				metapop_diploid_vector_type_allocator > * diploids,
  mutation_list_type<typename gamete_type::mutation_type,mutation_list_type_allocator > * mutations, 
	       const unsigned * N_curr, 
	       const unsigned * N_next, 
	       const double & mu,
	       const mutation_model & mmodel,
	       const recombination_policy & rec_pol,
	       const mutation_insertion_policy & mpolicy,
	       const gamete_insertion_policy & gpolicy_mut,
	       const diploid_fitness_function_container & ffs,
	       const mutation_removal_policy & mp,
	       const migration_policy & mig,
	       const double * f = NULL);

  /*! \brief recombination for individual-based forward simulations
    \param r GSL random number generator
    \param littler The probability of a single recombination event between g1 and g2
    \param gametes Pointer to the list of gametes in the population
    \param g1 Iterator to the first gamete involved in the recombination event
    \param g2 Iterator to the second gamete involved in the recombination event
    \param mf Recombination policy which generates crossover positions
    \param gpolicy Policy determining how new gametes are added to population

    \note g1 and g2 will be changed
    \note The type of g1 and g2 is gamete_list_type<gamete_type,list_type_allocator >::iterator
    \note The return value may be 0 even if littler is large.  The code recognizes when crossovers could not modify the gametes, and the function returns when such cases are found
    \return The number of crossovers that happened between g1 and g2 (which is Poisson with mean littler)
   */
  template< typename iterator_type,
	    typename gamete_insertion_policy,
	    typename recombination_map,
	    typename list_type_allocator,
	    template<typename,typename> class list_type>
  unsigned recombine_gametes( gsl_rng * r,
			      const double & littler,
			      list_type< typename iterator_type::value_type,list_type_allocator > * gametes,
			      iterator_type & g1,
			      iterator_type & g2,
			      const recombination_map & mf,
			      const gamete_insertion_policy & gpolicy);

  //Multilocus models
  //single deme, constant N
  template< typename gamete_type,
	    typename glist_vector_type_allocator,
	    typename gamete_list_type_allocator,
	    typename mutation_list_type_allocator,
	    typename diploid_vector_type_allocator,
	    typename locus_vector_type_allocator,
	    typename diploid_fitness_function,
	    typename mutation_removal_policy,
	    typename mutation_model_container,
	    typename recombination_policy_container,
	    typename mutation_insertion_policy,
	    typename gamete_insertion_policy,
	    template<typename,typename> class gamete_list_type,
	    template<typename,typename> class glist_vector_type,
	    template<typename,typename> class mutation_list_type,
	    template<typename,typename> class diploid_vector_type,
	    template<typename,typename> class locus_vector_type>
  double
  sample_diploid(gsl_rng * r,
		 glist_vector_type< gamete_list_type<gamete_type,
		 gamete_list_type_allocator> ,
		 glist_vector_type_allocator > * gametes,
		 diploid_vector_type<locus_vector_type<std::pair<typename gamete_list_type< gamete_type,gamete_list_type_allocator >::iterator,
		 typename gamete_list_type< gamete_type,gamete_list_type_allocator >::iterator>,
						       locus_vector_type_allocator>,
				     diploid_vector_type_allocator> * diploids,
		 mutation_list_type<typename gamete_type::mutation_type,mutation_list_type_allocator > * mutations, 
		 const unsigned & N
		 const double * mu,
		 const mutation_model_container & mmodel,
		 const recombination_policy_container & rec_policies,
		 const double * r_between_loci,
		 const mutation_insertion_policy & mpolicy,
		 const gamete_insertion_policy & gpolicy_mut,
		 const diploid_fitness_function & ff,
		 const mutation_removal_policy & mp,
		 const double & f);

//single deme, N changing
template< typename gamete_type,
	  typename glist_vector_type_allocator,
	  typename gamete_list_type_allocator,
	  typename mutation_list_type_allocator,
	  typename diploid_vector_type_allocator,
	  typename locus_vector_type_allocator,
	  typename diploid_fitness_function,
	  typename mutation_removal_policy,
	  typename mutation_model_container,
	  typename recombination_policy_container,
	  typename mutation_insertion_policy,
	  typename gamete_insertion_policy,
	  template<typename,typename> class gamete_list_type,
	  template<typename,typename> class glist_vector_type,
	  template<typename,typename> class mutation_list_type,
	  template<typename,typename> class diploid_vector_type,
	  template<typename,typename> class locus_vector_type>
double
sample_diploid(gsl_rng * r,
	       glist_vector_type< gamete_list_type<gamete_type,
						   gamete_list_type_allocator> ,
	       glist_vector_type_allocator > * gametes,
	       diploid_vector_type<locus_vector_type<std::pair<typename gamete_list_type< gamete_type,gamete_list_type_allocator >::iterator,
							       typename gamete_list_type< gamete_type,gamete_list_type_allocator >::iterator>,
						     locus_vector_type_allocator>,
	       diploid_vector_type_allocator> * diploids,
	       mutation_list_type<typename gamete_type::mutation_type,mutation_list_type_allocator > * mutations, 
	       const unsigned & N_curr, 
	       const unsigned & N_next, 
	       const double * mu,
	       const mutation_model_container & mmodel,
	       const recombination_policy_container & rec_policies,
	       const double * r_between_loci,
	       const mutation_insertion_policy & mpolicy,
	       const gamete_insertion_policy & gpolicy_mut,
	       const diploid_fitness_function & ff,
	       const mutation_removal_policy & mp,
	       const double & f);
}
#endif /* _DIPLOID_FUNCTIONS_HPP_ */
#include <fwdpp/diploid_functions.tcc>


