#ifndef _MUTATION_HPP_
#define _MUTATION_HPP_

#include <fwdpp/forward_types.hpp>
#include <fwdpp/fwd_functional.hpp>
#include <algorithm>
#include <limits>
#include <cmath>
#include <boost/bind.hpp>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

namespace KTfwd
{
  /*! \brief Apply mutation model to a population
    \example diploid.cc
    \param r GSL random number generator
    \param gametes Pointer to the vector of gametes in the population
    \param mutations Pointer to the list of mutations in the population
    \param mu the TOTAL mutation rate per gamete
    \param mmodel Mutation model policy
    \param gpolicy Policy determining how new gametes are added to population
    \param mpolicy Policy determining how new mutations are added to the population
   */
  template< typename gamete_type,
	    typename mutation_model,
	    typename gamete_insertion_policy,
	    typename mutation_insertion_policy,
	    typename vector_type_allocator,
	    typename list_type_allocator,
	    template<typename,typename> class vector_type,
	    template<typename,typename> class list_type>
  unsigned mutate(gsl_rng * r, 
		  vector_type<gamete_type,vector_type_allocator > * gametes, 
		  list_type<typename gamete_type::mutation_type,list_type_allocator > * mutations, 
		  const double & mu,
		  const mutation_model & mmodel,
		  const gamete_insertion_policy & gpolicy, 
		  const mutation_insertion_policy & mpolicy);

  /*! \brief Apply mutation model to an individual gamete.  Used for individual-based forward simulations
    \param r GSL random number generator
    \param gametes Pointer to the list of gametes in the population
    \param mutations Pointer to the list of mutations in the population
    \param g An iterator to the gamete that will be mutated by this function
    \param mu the TOTAL mutation rate per gamete
    \param mmodel Mutation model policy
    \param gpolicy Policy determining how new gametes are added to population
    \param mpolicy Policy determining how new mutations are added to the population

    \note g is passed non-const and will be modified by mutation events.
    \note The type of g is vector_type<gamete_type,vector_type_allocator >::iterator
    \note Used in invididual-based forward simulations.
    \return An iterator to the newly-created gamete, or to g if no mutation occurs.
   */
  template< typename iterator_type,
	    typename mutation_model,
	    typename mutation_insertion_policy,
	    typename gamete_insertion_policy,
	    typename list_type_allocator,
	    typename list_type_allocator2,
	    template<typename,typename> class list_type,
	    template<typename,typename> class list_type2>
  iterator_type mutate_gamete( gsl_rng * r,
			       const double & mu, list_type< typename iterator_type::value_type,list_type_allocator > * gametes,
			       list_type2<typename iterator_type::value_type::mutation_type,list_type_allocator2 > * mutations, 
			       iterator_type &g,
			       const mutation_model &mmodel,
			       const mutation_insertion_policy & mpolicy,
			       const gamete_insertion_policy & gpolicy);
}
#endif /* _MUTATION_HPP_ */
#include <fwdpp/mutation.tcc>


