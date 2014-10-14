#ifndef __FWDPP_GAMETE_BASED_HPP__
#define __FWDPP_GAMETE_BASED_HPP__

#include <vector>

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
}

#endif
