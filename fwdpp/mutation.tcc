//  -*- C++ -*- 
#ifndef _FWDPP_MUTATION_TCC_
#define _FWDPP_MUTATION_TCC_

#include <type_traits>
#include <algorithm>
#include <numeric>

#include <gsl/gsl_randist.h>

#include <fwdpp/internal/mutation_internal.hpp>
#include <vector>

namespace KTfwd
{
  template< typename queue_type,
	    typename queue_type2,
	    typename iterator_type,
	    typename mutation_model,
	    typename gamete_insertion_policy,
	    typename list_type_allocator,
	    typename list_type_allocator2,
	    template<typename,typename> class list_type,
	    template<typename,typename> class list_type2>
  iterator_type mutate_gamete_recycle( queue_type & recycling_bin,
				       queue_type2 & gamete_recycling_bin,
				       gsl_rng * r,
				       const double & mu, list_type< typename iterator_type::value_type,list_type_allocator > * gametes,
				       list_type2<typename iterator_type::value_type::mutation_type,list_type_allocator2 > * mutations, 
				       iterator_type &g,
				       const mutation_model &mmodel,
				       const gamete_insertion_policy & gpolicy)
  {
    assert( g != gametes->end() );
    unsigned nm = gsl_ran_poisson(r,mu);
    if ( nm )
      {
	assert( g->n > 0 );
	g->n--;
	//Recycle an already-allocated gamete, if possible
	if (!gamete_recycling_bin.empty())
	  {
	    auto __g = gamete_recycling_bin.front();
	    __g->mutations=g->mutations;
	    __g->smutations=g->smutations;
	    __g->n=1;
	    fwdpp_internal::add_N_mutations_recycle(recycling_bin,mmodel,nm,mutations,*__g);
	    gamete_recycling_bin.pop();
	    return __g;
	  }
	typename iterator_type::value_type ng( 1, g->mutations,g->smutations);
	fwdpp_internal::add_N_mutations_recycle(recycling_bin,mmodel,nm,mutations,ng);
	return gpolicy(std::move(ng),gametes);
      }
    return g;
  }
}
#endif /* _FWDPP_MUTATION_TCC_ */
