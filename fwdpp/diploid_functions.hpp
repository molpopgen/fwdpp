#ifndef _DIPLOID_FUNCTIONS_HPP_
#define _DIPLOID_FUNCTIONS_HPP_

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/*! \file diploid_functions.hpp
  \brief Wright-Fisher sampling in a finite population and recombination 
*/ 

//There are too many variants of sample_diploid to list here, so they are spawned off into other headers
#include <fwdpp/diploid_gamete_based.hpp>
#include <fwdpp/diploid_individual_based.hpp>
#include <fwdpp/diploid_individual_based_multilocus.hpp>

namespace KTfwd
{
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
	    typename recombination_map,
	    typename list_type_allocator,
	    template<typename,typename> class list_type>
  unsigned recombine_gametes( gsl_rng * r,
			      const double & littler,
			      list_type< typename iterator_type::value_type,list_type_allocator > * gametes,
			      iterator_type & g1,
			      iterator_type & g2,
			      const recombination_map & mf);


  //Multilocus models

}
#endif /* _DIPLOID_FUNCTIONS_HPP_ */
#include <fwdpp/diploid_functions.tcc>


