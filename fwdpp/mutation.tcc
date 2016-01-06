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
	    typename mutation_model,
	    typename gamete_insertion_policy,
	    typename gcont_t,
	    typename mcont_t>
  std::size_t mutate_gamete_recycle( queue_type & recycling_bin,
				     queue_type2 & gamete_recycling_bin,
				     gsl_rng * r,
				     const double & mu,
				     gcont_t & gametes,
				     mcont_t & mutations,
				     const size_t g,
				     const mutation_model &mmodel,
				     const gamete_insertion_policy & gpolicy)
  {
    assert(g<gametes.size());
    //assert( g != gametes->end() );
    unsigned nm = gsl_ran_poisson(r,mu);
    if ( nm )
      {
	assert(gametes[g].n);
	//assert( g->n > 0 );
	gametes[g].n--;
	//g->n--;
	//Recycle an already-allocated gamete, if possible
	if (!gamete_recycling_bin.empty())
	  {
	    auto __g = gamete_recycling_bin.front();
	    assert(!gametes[__g].n);
	    gametes[__g].mutations=gametes[g].mutations;
	    gametes[__g].smutations=gametes[g].smutations;
	    gametes[__g].n=1;
	    fwdpp_internal::add_N_mutations_recycle(recycling_bin,mmodel,nm,mutations,gametes[__g]);
	    assert( std::is_sorted( gametes[__g].mutations.begin(),
				    gametes[__g].mutations.end(),
				    [&mutations](const size_t i, const size_t j) {
				      return mutations[i].pos<mutations[j].pos;
				    } ) );
	    assert( std::is_sorted( gametes[__g].smutations.begin(),
				    gametes[__g].smutations.end(),
				    [&mutations](const size_t i, const size_t j) {
				      return mutations[i].pos<mutations[j].pos;
				    } ) );
	    gamete_recycling_bin.pop();
	    return __g;
	  }
	typename gcont_t::value_type ng( 1, gametes[g].mutations,gametes[g].smutations);
	fwdpp_internal::add_N_mutations_recycle(recycling_bin,mmodel,nm,mutations,ng);
	assert( std::is_sorted( ng.mutations.begin(),
				ng.mutations.end(),
				[&mutations](const size_t i, const size_t j) {
				  return mutations[i].pos<mutations[j].pos;
				} ) );
	assert( std::is_sorted( ng.smutations.begin(),
				ng.smutations.end(),
				[&mutations](const size_t i, const size_t j) {
				  return mutations[i].pos<mutations[j].pos;
				} ) );
	return gpolicy(std::move(ng),gametes);
      }
    return g;
  }
}
#endif /* _FWDPP_MUTATION_TCC_ */
