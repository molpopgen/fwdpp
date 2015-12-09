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
	     typename diploid_type_itr,
	     typename gamete_lookup_t,
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
	     typename bw_locus_rec_fxn
#endif
    >
    void multilocus_rec_mut(gsl_rng * r,
			    diploid_type parent1, //Copy--this is intentional
			    diploid_type parent2, //Copy--this is intentional
			    diploid_type_itr & offspring, //non-const ref, again intentional
			    mqueue_t & mutation_recycling_bin,
			    gqueue_t & gamete_recycling_bin,
			    gamete_lookup_t & gamete_lookup,
			    const recombination_policy_container & rec_pols,
			    const bw_locus_rec_fxn & blrf,
			    const double * r_bw_loci,
			    const int iswitch1,
#ifndef FWDPP_UNIT_TESTING
			    const int iswitch2,
			    glist_t * gametes,
			    mlist_t * mutations,
			    const double * mu,
			    const mutation_model_container & mmodel,
			    const gamete_insertion_policy & gpolicy_mut
#else
			    const int iswitch2
#endif			
			    )
    {
      //I see the problem: how to get the positions ahead of time...
      //Maybe we can simply increment all downstream values by 1 if a swap is needed, and do so if odd?

      std::vector<int> nswaps1(parent1.size(),iswitch1),nswaps2(parent2.size(),iswitch2);
      std::vector<int>::iterator s1 = nswaps1.begin(),s2=nswaps2.begin();
      typename diploid_type_itr::value_type::iterator optr = offspring->begin();
      for( unsigned i = 0 ; i < parent1.size() ; ++i,++s1,++s2,++optr )
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

	  //Assign pointers to offspring
	  optr->first = parent1[i].first;
	  optr->second = parent2[i].first;

	  //mechanics of recombination
	  unsigned nrec = rec_pols[i](optr->first,parent1[i].second,gamete_lookup,gamete_recycling_bin);
	  if(nrec%2!=0.) std::transform( s1+1,nswaps1.end(),s1+1,std::bind(std::plus<int>(),std::placeholders::_1,nrec) );

	  nrec = rec_pols[i](optr->second,parent2[i].second,gamete_lookup,gamete_recycling_bin);
	  if(nrec%2!=0.) std::transform( s2+1,nswaps2.end(),s2+1,std::bind(std::plus<int>(),std::placeholders::_1,nrec) );

#ifndef FWDPP_UNIT_TESTING
	  optr->first->n++;
	  optr->second->n++;
	  optr->first = mutate_gamete_recycle(mutation_recycling_bin,gamete_recycling_bin,r,mu[i],gametes,mutations,optr->first,mmodel[i],gpolicy_mut);
	  optr->second = mutate_gamete_recycle(mutation_recycling_bin,gamete_recycling_bin,r,mu[i],gametes,mutations,optr->second,mmodel[i],gpolicy_mut);
#endif
	}
    }
  }
}

#endif
