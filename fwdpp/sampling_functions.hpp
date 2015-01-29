#ifndef __SAMPLING_FUNCTIONS_HPP__
#define __SAMPLING_FUNCTIONS_HPP__

#include <vector>
#include <map>
#include <functional>
#include <algorithm>
#include <utility>
#include <string>
#include <cassert>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/*! @defgroup samplingPops Functions related to taking samples from simulated populations

  This collection of functions allows a user to draw a sample of size \f$n \ll 2N\f$ from
  a simulated population.

  The library provides several overloads of the functions KTfwd::ms_sample and KTfwd::ms_sample_separate.  Further,
  the behavior of these functions differs slightly dpeending on whether or not the output from the simulation 
  comes from a gamete- or an indvidual-based method.
*/

/*!
  @defgroup samplingPopsGamete Gamete-based simulation output
  @ingroup samplingPops
*/

/*!
  @defgroup samplingPopsInd Individual-based simulation output
  @ingroup samplingPops
*/

namespace KTfwd
{
  /* \brief Site frequency spectrum
     \return Site frequency spectrum
  */
  template<typename iterator_type>
  std::vector<unsigned> population_sfs( iterator_type beg,
					iterator_type end,
					const unsigned & N)
  {
    std::vector<unsigned> psfs(N-1,0);
    while(beg != end)
      {
	if(beg->n >0 && beg->n < N) psfs[beg->n-1]++;
	beg++;
      }
    return psfs;
  }


  template< typename gamete_type,
	    typename vector_type_allocator,
	    template<typename,typename> class vector_type>
  std::vector<unsigned> sample_sfs(gsl_rng * r, 
				   const vector_type<gamete_type,vector_type_allocator > & gametes,
				   const unsigned & n, const unsigned & N);
  /*!
    Take a sample of size n from a larger population of N diploids
    \ingroup samplingPopsGamete
  */
  template< typename gamete_type,
	    typename vector_type_allocator,
	    template<typename,typename> class vector_type>
  std::vector< std::pair<double, std::string> > 
  ms_sample(gsl_rng * r,
	    const vector_type<gamete_type,vector_type_allocator > * gametes,
	    const unsigned & n, const unsigned & N,
	    bool remove_fixed = true);
  
  
  /*!
    \ingroup samplingPopsGamete
  */
  template< typename gamete_type,
	    typename vector_type_allocator,
	    template<typename,typename> class vector_type>
  std::pair< std::vector< std::pair<double, std::string> > ,
	     std::vector< std::pair<double, std::string> > >
  ms_sample_separate(gsl_rng * r,
		     const vector_type<gamete_type,vector_type_allocator > * gametes,
		     const unsigned & n, const unsigned & N,
		     bool remove_fixed = true);

  /*!
    \brief Sampling from a population in an individual-based simulation
    \ingroup samplingPopsInd
  */
  template<typename iterator_type,
	   typename allocator,
	   template<typename,typename> class vector_type >
  std::vector< std::pair<double,std::string> >
  ms_sample( gsl_rng * r,
	     const vector_type< std::pair<iterator_type,iterator_type>, allocator > * diploids,
	     const unsigned & n,
	     const bool & remove_fixed = true);

  /*!
    \brief Sampling from a population in an individual-based simulation.  Selected and neutral mutations returned separately
    \ingroup samplingPopsInd
  */
  template<typename iterator_type,
	   typename allocator,
	   template<typename,typename> class vector_type >
  std::pair<std::vector< std::pair<double,std::string> >,
	    std::vector< std::pair<double,std::string> > >
  ms_sample_separate( gsl_rng * r,
		      const vector_type< std::pair<iterator_type,iterator_type>, allocator > * diploids,
		      const unsigned & n,
		      const bool & remove_fixed = true);

  /*!
    \brief Sample from an individual-based, multi-locus simulation.
    \return A vector of vectors of variable sites.  There is 1 vector per locus.
    \note Neutral + selected mutations intermixed
    \ingroup samplingPopsInd
  */
  template<typename iterator_type,
	   typename allocator,
	   typename outer_allocator,
	   template<typename,typename> class vector_type,
	   template<typename,typename> class outer_vector_type>
  std::vector< std::vector< std::pair<double,std::string> > >
  ms_sample( gsl_rng * r,
	     const outer_vector_type< vector_type< std::pair<iterator_type,iterator_type>, allocator >, outer_allocator > * diploids,
	     const unsigned & n,
	     const bool & remove_fixed);

  /*!
    \brief Sample from an individual-based, multi-locus simulation.
    \return A vector of pairs of vectors of variable sites.  There is 1 vector per locus.
    \note For each locus, the first member of the pair corresponds to neutral sites, the second to selected.
    \ingroup samplingPopsInd
  */
  template<typename iterator_type,
	   typename allocator,
	   typename outer_allocator,
	   template<typename,typename> class vector_type,
	   template<typename,typename> class outer_vector_type>
  std::vector< std::pair<std::vector< std::pair<double,std::string> >,std::vector< std::pair<double,std::string> > > >
  ms_sample_separate( gsl_rng * r,
		      const outer_vector_type< vector_type< std::pair<iterator_type,iterator_type>, allocator >, outer_allocator > * diploids,
		      const unsigned & n,
		      const bool & remove_fixed);
}
#endif 
#include <fwdpp/sampling_functions.tcc>
