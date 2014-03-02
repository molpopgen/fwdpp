#ifndef __RECOMBINATION_METHODS_IND_HPP__
#define __RECOMBINATION_METHODS_IND_HPP__

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

namespace KTfwd
{

struct genetics101
/*! Genetics 101: simple model of recombination.  r is the probability that the two gametes recombine
*/
{
  typedef unsigned result_type;
  template<typename gamete_iterator_type,
	   typename gamete_list_type_allocator,
	   template<typename,typename> class gamete_list_type,
	   typename rec_pos_generator>
  unsigned operator()( gamete_iterator_type & g1,
		       gamete_iterator_type & g2,
		       gamete_list_type< typename gamete_iterator_type::value_type, gamete_list_type_allocator > * gametes,
		       const double & littler,
		       gsl_rng * r,
		       const rec_pos_generator & rp) const
  {
    typedef gamete_list_type< typename gamete_iterator_type::value_type, gamete_list_type_allocator > glist_t;
    unsigned NREC = 0;
    if( g1 != g2 )
      //then a non-parental type is inherited from p1 and p1 has two different gametes
      {
	NREC += recombine_gametes(r,littler,gametes,g1,g2,rp,
				  boost::bind(update_if_exists_insert<typename gamete_iterator_type::value_type,glist_t>,_1,gametes));	  
      }
    return NREC;
  }	   
};

}

#endif //__RECOMBINATION_METHODS_IND_HPP__
