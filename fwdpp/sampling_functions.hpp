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
  the behavior of these functions differs slightly depending on whether or not the output from the simulation 
  comes from a gamete- or an indvidual-based method.

  The following features are in common to all versions of these functions:

  1. They only support bi-allelic mutation positions.
  2. All functions return some number of objects of type std::vector< std::pair<double, std::string> >.  
  Each element in the vector corresponds to a mutation.  The double corresponds to the mutation position,
  and the string corresponds to the character states in each of the \f$n\f$ samples.  The strings contain the character
  '0' when a sample has the ancestral state or '1' for the derived (mutant) state.  This type may be used
  to populate a polymorphism table from the [libsequence](http://molpopgen.github.io/libsequence/) library:
  \code
  auto x = KTfwd::ms_sample( appropriate arguments );
  Sequence::SimData xx(x.begin(),x.end());
  \endcode
  When a SimData object is written to a stream, its output format is the same as that used by Dick Hudson's coalescent simulation
  program [ms](http://home.uchicago.edu/~rhudson1/source/mksamples.html)
  3. All versions of KTfwd::ms_sample return vectors where the mutations affecting fitness and those not affecting fitness are intermingled.
  4. All versions of KTfwd::ms_sample_separate return pairs of vectors separating the mutations not affecting fitness from those that do.  
  The first member of each pair is a vector of "neutral" mutations, and the second member is the vector of "selected" mutations.
  5. The functions are all implemented using sampling with replacement from the population.  Thus, setting \f$n = N\f$ (or \f$2N\f$) will NOT
  return all the variable sites in the entire population!!!

  A vector may contain data looking like the following:
  \verbatim
  0.0760 "10000"
  0.6805 "00010"
  \endverbatim

  The two rows are the two different positions (0.076 and 0.6805).  The first haplotype is "10", the second is "00", and the fourth is "01".
  
  If you converted the data to a Sequence::SimData object and printed it to screen or a file, the output would be:

  \verbatim
  //
  segsites: 2
  positions: 0.0760 0.6805 
  10
  00
  00
  01
  00
  \endverbatim
*/

/*!
  @defgroup samplingPopsGamete Randomly-sampling individual gametes.
  @ingroup samplingPops

  These functions draws a sample of size \f$n \ll 2N\f$ from a simulated population.  

  The gametes are randomly-sampled (with replacement) proportionally to their frequency in the population.

  The object passed to these functions is the container of gametes (e.g, some container of type KTfwd::gamete_base).

  The return values have no relation to any actual diploid individual in the population. Each haplotype
  corresponds to a simulated haplotype, but adjacent haplotypes are basically random draws of M&Ms from a jar.

  For these functions, \f$n\f$ can be odd or even.
*/

/*!
  @defgroup samplingPopsInd Randomly-sampling diploids.
  @ingroup samplingPops

  These functions pull a sample of \f$n\f$ simulated individuals from the population.  You may only use them on the 
  output of individual-based simulations (gamete-based simulations never explicitly store individuals).

  The individuals are sampled uniformly, and with replacement, with no regard for their fitness.

  The object passed to these functions is a container of diploids. 

  For these functions, \f$n\f$ must be an even number, as it represents the number of alleles to sample (twice the 
  number of individuals).

  The return values of these functions store individuals in the order that they were sampled.
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
	    typename allocator_t,
	    template<typename,typename> class container_t>
  std::vector<unsigned> sample_sfs(gsl_rng * r, 
				   const container_t<gamete_type,allocator_t > & gametes,
				   const unsigned & n, const unsigned & N);
  /*!
    Take a sample of size n from a larger population of N diploids
    \ingroup samplingPopsGamete
  */
  template< typename gamete_type,
	    typename allocator_t,
	    template<typename,typename> class container_t>
  //std::vector< std::pair<double, std::string> >
  typename std::enable_if< std::is_base_of<mutation_base,typename gamete_type::mutation_type>::value,
			   std::vector< std::pair<double, std::string> > >::type
  ms_sample(gsl_rng * r,
	    const container_t<gamete_type,allocator_t > * gametes,
	    const unsigned & n, const unsigned & N,
	    bool remove_fixed = true);
  
  
  /*!
    \ingroup samplingPopsGamete
  */
  template< typename gamete_type,
	    typename allocator_t,
	    template<typename,typename> class container_t>
  //std::pair< std::vector< std::pair<double, std::string> > ,
  //	     std::vector< std::pair<double, std::string> > >
  typename std::enable_if< std::is_base_of<mutation_base,typename gamete_type::mutation_type>::value,
			   std::pair< std::vector< std::pair<double, std::string> > ,
				      std::vector< std::pair<double, std::string> > > >::type
  ms_sample_separate(gsl_rng * r,
		     const container_t<gamete_type,allocator_t > * gametes,
		     const unsigned & n, const unsigned & N,
		     bool remove_fixed = true);

  /*!
    \brief Sampling from a population in an individual-based simulation
    \ingroup samplingPopsInd
  */
  template<//typename iterator_type,
	   typename allocator,
	   typename diploid_t,
	   template<typename,typename> class vector_type >
  //std::vector< std::pair<double,std::string> >
  typename std::enable_if< std::is_base_of<mutation_base,typename diploid_t::first_type::value_type::mutation_type>::value,
			   std::vector< std::pair<double,std::string> > >::type
  ms_sample( gsl_rng * r,
	     const vector_type< diploid_t, allocator > * diploids,
	     const unsigned & n,
	     const bool & remove_fixed = true);

  /*!
    \brief Sampling from a population in an individual-based simulation.  Selected and neutral mutations returned separately
    \ingroup samplingPopsInd
  */
  template<//typename iterator_type,
	   typename allocator,
	   typename diploid_t,
	   template<typename,typename> class vector_type >
  std::pair<std::vector< std::pair<double,std::string> >,
	    std::vector< std::pair<double,std::string> > >
  ms_sample_separate( gsl_rng * r,
		      const vector_type< diploid_t, allocator > * diploids,
		      const unsigned & n,
		      const bool & remove_fixed = true);

  /*!
    \brief Sample from an individual-based, multi-locus simulation.
    \return A vector of vectors of variable sites.  There is 1 vector per locus.
    \note Neutral + selected mutations intermixed
    \ingroup samplingPopsInd
  */
  template<typename diploid_type,
	   typename allocator,
	   typename outer_allocator,
	   template<typename,typename> class vector_type,
	   template<typename,typename> class outer_vector_type>
  std::vector< std::vector< std::pair<double,std::string> > >
  ms_sample( gsl_rng * r,
	     const outer_vector_type< vector_type< diploid_type, allocator >, outer_allocator > * diploids,
	     const unsigned & n,
	     const bool & remove_fixed);

  /*!
    \brief Sample from an individual-based, multi-locus simulation.
    \return A vector of pairs of vectors of variable sites.  There is 1 vector per locus.
    \note For each locus, the first member of the pair corresponds to neutral sites, the second to selected.
    \ingroup samplingPopsInd
  */
  template<typename diploid_type,
	   typename allocator,
	   typename outer_allocator,
	   template<typename,typename> class vector_type,
	   template<typename,typename> class outer_vector_type>
  //std::vector< std::pair<std::vector< std::pair<double,std::string> >,std::vector< std::pair<double,std::string> > > >
  typename std::enable_if< std::is_base_of<mutation_base,typename diploid_type::first_type::value_type::mutation_type>::value,
			   std::pair<std::vector< std::pair<double,std::string> >,
				     std::vector< std::pair<double,std::string> > > >::type
  ms_sample_separate( gsl_rng * r,
		      const outer_vector_type< vector_type< diploid_type, allocator >, outer_allocator > * diploids,
		      const unsigned & n,
		      const bool & remove_fixed);
}
#endif 
#include <fwdpp/sampling_functions.tcc>
