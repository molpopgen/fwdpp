#ifndef __FWDPP_GENETICS101_HPP__
#define __FWDPP_GENETICS101_HPP__

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

namespace KTfwd
{
  struct genetics101
  /*! Genetics 101: simple model of recombination.  r is the probability that the two gametes recombine
   */
  {
    using result_type = unsigned;
    template<typename gamete_iterator_type,
	     typename gamete_list_type_allocator,
	     typename glookup_t,
	     typename queue_t,
	     template<typename,typename> class gamete_list_type,
	     typename rec_pos_generator>
    unsigned operator()( gamete_iterator_type & g1,
			 gamete_iterator_type & g2,
			 glookup_t & gamete_lookup,
			 queue_t & gamete_recycling_bin,
			 typename gamete_iterator_type::value_type::mutation_container & neutral,
			 typename gamete_iterator_type::value_type::mutation_container & selected,
			 gamete_list_type< typename gamete_iterator_type::value_type, gamete_list_type_allocator > * gametes,
			 const double & littler,
			 gsl_rng * r,
			 const rec_pos_generator & rp) const
    {
      unsigned NREC = 0;
      if( g1 != g2 )
	//then a non-parental type is inherited from p1 and p1 has two different gametes
	{
	  NREC += recombine_gametes(r,littler,gametes,g1,g2,gamete_lookup,gamete_recycling_bin,neutral,selected,rp);
	}
      return NREC;
    }
  };
}

#endif //__FWDPP_GENETICS101_HPP__
