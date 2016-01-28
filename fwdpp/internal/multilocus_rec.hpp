#ifndef __FWDPP_INTERNAL_MULTILOCUS_REC_HPP__
#define __FWDPP_INTERNAL_MULTILOCUS_REC_HPP__

/*
  The mechanics of crossing over for a multilocus
  simulation
 */
#include <vector>
#include <algorithm>
#include <gsl/gsl_rng.h>

namespace KTfwd {
  namespace fwdpp_internal {

    /*!
      Mechanics of segregation, recombination, and mutation for multi-locus API
      This template declaratiopn is messy due to the FWDPP_UNIT_TESTING symbol.

      The reason is that I want to be able to unit test the recombination portion
      of this code in isolation, w/o having to do anything involving mutation policies.

      Thus, the relevant unit tests define FWDPP_UNIT_TESTING so that the mutation-related
      stuff is skipped during compilation.
    */
    template<typename diploid_type,
	     typename recombination_policy_container,
	     typename mqueue_t,
	     typename gqueue_t,
#ifndef FWDPP_UNIT_TESTING
	     typename bw_locus_rec_fxn,
	     typename mlist_t,
	     typename glist_t,
	     typename mutation_model_container,
	     typename gamete_insertion_policy
#else
	     typename mlist_t,
	     typename glist_t,
	     typename bw_locus_rec_fxn
#endif
    >
    diploid_type multilocus_rec_mut(gsl_rng * r,
				    diploid_type parent1,
				    diploid_type parent2,
				    mqueue_t & mutation_recycling_bin,
				    gqueue_t & gamete_recycling_bin,
				    const recombination_policy_container & rec_pols,
				    const bw_locus_rec_fxn & blrf,
				    const double * r_bw_loci,
				    const int iswitch1,
#ifndef FWDPP_UNIT_TESTING
				    const int iswitch2,
				    glist_t & gametes,
				    mlist_t & mutations,
				    typename glist_t::value_type::mutation_container & neutral,
				    typename glist_t::value_type::mutation_container & selected,
				    const double * mu,
				    const mutation_model_container & mmodel,
				    const gamete_insertion_policy & gpolicy_mut
#else
				    const int iswitch2,
				    glist_t & gametes,
				    mlist_t & mutations,
				    typename glist_t::value_type::mutation_container & neutral,
				    typename glist_t::value_type::mutation_container & selected
#endif
				    )
    {
      //I see the problem: how to get the positions ahead of time...
      //Maybe we can simply increment all downstream values by 1 if a swap is needed, and do so if odd?
      diploid_type offspring(parent1.size());
      std::vector<int> nswaps1(parent1.size(),iswitch1),nswaps2(parent2.size(),iswitch2);
      std::vector<int>::iterator s1 = nswaps1.begin(),s2=nswaps2.begin();

      for(unsigned i = 0 ; i < parent1.size() ; ++i,++s1,++s2)
	{
	  if(i)
	    {
	      // between-locus rec, parent 1
	      unsigned nrbw = blrf(r,r_bw_loci[i-1]);
	      //only modify if odd
	      if(nrbw%2!=0.) std::transform( s1,nswaps1.end(),s1,std::bind(std::plus<int>(),std::placeholders::_1,nrbw) );

	      // between-locus rec, parent 2
	      nrbw = blrf(r,r_bw_loci[i-1]);
	      //only modify if odd
	      if(nrbw%2!=0.) std::transform( s2,nswaps2.end(),s2,std::bind(std::plus<int>(),std::placeholders::_1,nrbw) );
	    }
	  //if ttl # recs before now is odd, swap parental pointers
	  if( *s1 % 2 != 0.) std::swap(parent1[i].first,parent1[i].second);
	  if( *s2 % 2 != 0.) std::swap(parent2[i].first,parent2[i].second);

	  //mechanics of recombination
	  auto xx = recombination(gametes,gamete_recycling_bin,neutral,selected,rec_pols[i],
				  parent1[i].first,parent1[i].second,mutations);
	  offspring[i].first = xx.first;
	  if(xx.second%2!=0.) std::transform( s1+1,nswaps1.end(),s1+1,std::bind(std::plus<int>(),std::placeholders::_1,xx.second) );

	  xx = recombination(gametes,gamete_recycling_bin,neutral,selected,rec_pols[i],
			    parent2[i].first,parent2[i].second,mutations);
	  offspring[i].second = xx.first;
	  if(xx.second%2!=0.) std::transform( s2+1,nswaps2.end(),s2+1,std::bind(std::plus<int>(),std::placeholders::_1,xx.second) );

#ifndef FWDPP_UNIT_TESTING
	  gametes[offspring[i].first].n++;
	  gametes[offspring[i].second].n++;
	  offspring[i].first = mutate_gamete_recycle(mutation_recycling_bin,gamete_recycling_bin,r,mu[i],gametes,mutations,offspring[i].first,mmodel[i],gpolicy_mut);
	  offspring[i].second = mutate_gamete_recycle(mutation_recycling_bin,gamete_recycling_bin,r,mu[i],gametes,mutations,offspring[i].second,mmodel[i],gpolicy_mut);
#endif
	}
      return offspring;
    }
  }
}

#endif
