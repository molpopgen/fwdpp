#ifndef __FWDPP_RECOMBINATION_HPP__
#define __FWDPP_RECOMBINATION_HPP__

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

namespace KTfwd
{
  struct poisson_xover
  {
    template<typename gamete_type,
	     typename mcont_t>
    std::vector<double> operator()(gsl_rng * r,
				   const double littler,
				   const double minpos,
				   const double maxpos,
				   //args below will be passed in by sample_diploid,
				   //but they are not used/needed here
				   //Replace
				   const gamete_type & ,
				   const gamete_type & ,
				   const mcont_t & ) const
    {
      unsigned nbreaks = (littler > 0) ? gsl_ran_poisson(r,littler) : 0u;
      if(!nbreaks) return {};

      std::vector<double> pos;
      pos.reserve(nbreaks+1);
      for(unsigned i = 0 ; i < nbreaks ; ++i)
	{
	  pos.emplace_back(gsl_ran_flat(r,minpos,maxpos));
	}
      std::sort(pos.begin(),pos.end());
      pos.emplace_back(std::numeric_limits<double>::max());
      return pos;
    }
  };

  template<typename gcont_t,
	   typename mcont_t,
	   typename lookup_t,
	   typename recbin_t,
	   typename recpol_t>
  std::size_t recombination(gcont_t & gametes,
			    lookup_t & gamete_lookup,
			    recbin_t & gamete_recycling_bin,
			    typename gcont_t::value_type::mutation_container & neutral,
			    typename gcont_t::value_type::mutation_container & selected,
			    const recpol_t & rec_pol,
			    const std::size_t g1,
			    const std::size_t g2,
			    const mcont_t & mutations);
    
  /*!
    Overload for fixed xover positions.
    Typically, this is called by the version taking 
    a recombination policy as an argument.

    If you wish to call this version directly,
    only do so if length pos > 1, pos is sorted
    in ascending order, and the last value in 
    pos is std::numeric_limits<double>::max(),
    which is assumed to be a terminating value 
    larger than any possible value for a mutation's position.

    \param pos A vector (with interface of std::vector) containing recombination breakpoints.  See note below.
    \param gametes A container of the gametes segregating in the population
    \param g1 An iterator, derived from gametes, representing one parental gamete.
    \param g2 An iterator, derived from gametes, representing the other parental gamete.
    \param gamete_lookup The return value of KTfwd::fwdpp_internal::gamete_lookup_table, which gets passed in via sample_diploid
    \return The number of breakpoints, which equals pos.size() - 1, as that is fixed in this case.
    \note The vector pos must be sorted (ascending order) and must contain the value std::numeric_limits<double>::max() as a terminating value.
  */
  template< typename iterator_type,
  	    typename vec_t,
  	    typename glist_t,
  	    typename glookup_t,
  	    typename queue_t>
  unsigned recombine_gametes( const vec_t & pos,
  			      glist_t * gametes,
  			      iterator_type & g1,
  			      iterator_type & g2,
  			      glookup_t & gamete_lookup,
  			      queue_t & gamete_recycling_bin,
  			      typename iterator_type::value_type::mutation_container & neutral,
  			      typename iterator_type::value_type::mutation_container & selected );

}
#endif // __FWDPP_RECOMBINATION_HPP__ 
#include <fwdpp/recombination.tcc>


