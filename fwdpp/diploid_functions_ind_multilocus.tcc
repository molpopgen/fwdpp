//  -*- C++ -*- 
#ifndef __DIPLOID_FUNCTIONS_IND_MULTILOC_TCC__
#define __DIPLOID_FUNCTIONS_IND_MULTILOC_TCC__

#include <fwdpp/mutation.hpp>
#include <fwdpp/internal/gsl_discrete.hpp>
#include <fwdpp/internal/multilocus_rec.hpp>

namespace KTfwd
{
  //single deme, N changing
  template< typename gamete_type,
	    typename glist_vector_type_allocator,
	    typename gamete_list_type_allocator,
	    typename mutation_list_type_allocator,
	    typename diploid_vector_type_allocator,
	    typename locus_vector_type_allocator,
	    typename diploid_fitness_function,
	    typename mutation_removal_policy,
	    typename mutation_model_container,
	    typename recombination_policy_container,
	    typename mutation_insertion_policy,
	    typename gamete_insertion_policy,
	    typename bw_locus_rec_fxn,
	    template<typename,typename> class gamete_list_type,
	    template<typename,typename> class glist_vector_type,
	    template<typename,typename> class mutation_list_type,
	    template<typename,typename> class diploid_vector_type,
	    template<typename,typename> class locus_vector_type>
  double
  sample_diploid(gsl_rng * r,
		 glist_vector_type< gamete_list_type<gamete_type,
		 gamete_list_type_allocator> ,
		 glist_vector_type_allocator > * gametes,
		 diploid_vector_type<locus_vector_type<std::pair<typename gamete_list_type< gamete_type,gamete_list_type_allocator >::iterator,
		 typename gamete_list_type< gamete_type,gamete_list_type_allocator >::iterator>,
		 locus_vector_type_allocator>,
		 diploid_vector_type_allocator> * diploids,
		 mutation_list_type<typename gamete_type::mutation_type,mutation_list_type_allocator > * mutations, 
		 const unsigned & N_curr, 
		 const unsigned & N_next, 
		 const double * mu,
		 const mutation_model_container & mmodel,
		 const recombination_policy_container & rec_policies,
		 const double * r_between_loci,
		 const bw_locus_rec_fxn & blrf,
		 const mutation_insertion_policy & mpolicy,
		 const gamete_insertion_policy & gpolicy_mut,
		 const diploid_fitness_function & ff,
		 const mutation_removal_policy & mp,
		 const double & f)
  {
#ifndef NDEBUG
    assert(diploids->size()==N_curr);
#endif

    //Adjust all mutation frequencies down to 0
    for( auto itr = mutations->begin() ;
	 itr != mutations->end() ; ++itr )
      {
	itr->n = 0;
      }
  
    //Vector of parental fitnesses
    std::vector<double> fitnesses(N_curr);
    double wbar = 0.;
  
    //Go over parents
    auto dptr = diploids->begin();
    for( unsigned i = 0 ; i < N_curr ; ++i,++dptr )
      {
	//set parental gamete counts to 0 for each locus
	for( auto j = dptr->begin() ; j != dptr->end() ; ++j )
	  {
	    j->first->n=0;
	    j->second->n=0;
	  }
	//Calculate the fitness of this parent
	fitnesses[i] += ff( *dptr );
	//increment pop. mean fitness
	wbar += fitnesses[i];
      }

#ifndef NDEBUG
    /*
      If we are debugging, let's make sure that every gamete has been set to n = 0.
      Rationale for check:  if we are failing to update data types properly, then
      it is possible that the "gamete pool" contains items not carried by any diploids.
      If so, this assertion will fail.
    */
    for ( auto i = gametes->cbegin() ; i != gametes->cend() ; ++i )
      {
	for (auto gptr = i->cbegin() ; gptr != i->cend() ; ++gptr )
	  {
	    assert( ! gptr->n );
	  }
      }
#endif

    wbar /= double(diploids->size());
    dptr = diploids->begin();  //reset to beginning of diploids
  
    fwdpp_internal::gsl_ran_discrete_t_ptr lookup(gsl_ran_discrete_preproc(fitnesses.size(),&fitnesses[0]));
 
    auto parents(*diploids); //Copy the parents.  Exact copy of diploids--same fitnesses, etc.
    const auto pptr = parents.cbegin();
 
    //Change the population size and reset dptr to avoid iterator invalidation
    if( diploids->size() != N_next)
      {
	diploids->resize(N_next);
	dptr = diploids->begin();
      }
  
    assert(diploids->size()==N_next);
 
    for( unsigned curr_dip = 0 ; curr_dip < N_next ; ++curr_dip )
      {
	assert(dptr==diploids->begin());
	assert( (dptr+curr_dip) < diploids->end() );

	//Choose the two parents
	typename decltype(pptr)::difference_type p1 = decltype(p1)(gsl_ran_discrete(r,lookup.get()));
	decltype(p1) p2  = (gsl_rng_uniform(r) <= f) ? p1 : decltype(p1)(gsl_ran_discrete(r,lookup.get()));
	assert(p1<N_curr);
	assert(p2<N_curr);
     
	//Need to store a vector of the equivalent of p1g1,p1g2 out to p1gn,p2gn
	//This is a trivial copying of iterators, so not that expensive
	auto p1c = *(pptr+p1),
	  p2c( *(pptr+p2) );
	/*
	  Through 0.2.9, we just said assert(p1c == *(pptr+p1)) here.
	  Oddly, that assert always fails if the type of a diploid is
	  boost::container::vector< std::pairs of iterators to gametes >.
	*/
#ifndef NDEBUG
	for(unsigned i = 0 ; i < p1c.size() ; ++i )
	  {
	    assert( (p1c.begin()+i)->first == ((pptr+p1)->begin()+i)->first );
	    assert( (p1c.begin()+i)->second == ((pptr+p1)->begin()+i)->second );
	  }
	for(unsigned i = 0 ; i < p2c.size() ; ++i )
	  {
	    assert( (p2c.begin()+i)->first == ((pptr+p2)->begin()+i)->first );
	    assert( (p2c.begin()+i)->second == ((pptr+p2)->begin()+i)->second );
	  }
#endif

	//Using just the routines below give correct E[S] for n = 20.
	//recombination, too! OK--this seems good based on some limited testing.
	//All testing so far based on n=2 using strobeck_morgan.cc

	//For 3 locus (test/strobeck_morgan.cc), the dist. of 
	//summary stats looks great vis-a-vis ms for N=20,theta=20,rho=1, N=1e3,10N gens
	bool p1g1 = (gsl_rng_uniform(r) <= 0.5) ? true : false,
	  p2g1 = (gsl_rng_uniform(r) <= 0.5) ? true : false;
	auto ptr2cdip = (dptr+curr_dip)->begin();
	bool LO1 = true, LO2 = true;
	for ( unsigned i = 0 ; i < p1c.size() ; ++i )
	  {
	    //This entire bit from here...
	    (ptr2cdip+i)->first = fwdpp_internal::multilocus_rec( r,rec_policies[i],blrf,
								  r_between_loci,i,
								  p1c[i].first,p1c[i].second,
								  p1g1,LO1 );
	    (ptr2cdip+i)->second = fwdpp_internal::multilocus_rec( r,rec_policies[i],blrf,
								   r_between_loci,i,
								   p2c[i].first,p2c[i].second,
								   p2g1,LO2 );
	    (ptr2cdip+i)->first->n++;
	    (ptr2cdip+i)->second->n++;

	    (ptr2cdip+i)->first = mutate_gamete( r,mu[i],&*(gametes->begin()+i),mutations,(ptr2cdip+i)->first,mmodel[i],mpolicy,gpolicy_mut);
	    (ptr2cdip+i)->second = mutate_gamete( r,mu[i],&*(gametes->begin()+i),mutations,(ptr2cdip+i)->second,mmodel[i],mpolicy,gpolicy_mut);	 
	  }
      }

    //update mutation counts in gametes
    auto glist_updater = []( decltype( *(gametes->begin()) ) & __g) {
      auto __first=__g.begin(),__last=__g.end();
      decltype(__first) __temp;
      while(__first!=__last)
	{
	  if(! __first->n)
	    {
	      __temp=__first;
	      ++__first;
	      __g.erase(__temp);
	    }
	  else
	    {
	      adjust_mutation_counts(__first,__first->n);
	      ++__first;
	    }
	}
    };
    std::for_each( gametes->begin(), gametes->end(), std::cref(glist_updater) );
    std::for_each(gametes->begin(),gametes->end(),
		  [&mp](decltype( *(gametes->begin()) ) & g) {
		    std::for_each(g.begin(),g.end(),
				  [&mp](decltype( *(g.begin()) ) & __g) {
				    __g.mutations.erase( std::remove_if(__g.mutations.begin(),__g.mutations.end(),std::cref(mp)),__g.mutations.end());
				    __g.smutations.erase( std::remove_if(__g.smutations.begin(),__g.smutations.end(),std::cref(mp)),__g.smutations.end());
				  });
		  });
#ifndef NDEBUG
    for ( auto i = gametes->begin() ; i != gametes->end() ; ++i )
      {
	unsigned sum = 0;
	for (auto gptr = i->begin() ; gptr != i->end() ; ++gptr )
	  {
	    //make sure that mutation frequencies are >= gamete frequencies
	    for( decltype(gptr->mutations.size()) j = 0 ; j < gptr->mutations.size() ; ++j )
	      {
		assert( gptr->mutations[j]->n >= gptr->n );
	      }
	    for( decltype(gptr->smutations.size()) j = 0 ; j < gptr->smutations.size() ; ++j )
	      {
		assert( gptr->smutations[j]->n >= gptr->n );
	      }
	    
	    sum += gptr->n;
	    assert( sum && sum <= 2*N_next );
	  }
	assert(sum == 2*N_next);
      }
#endif
    return wbar;
  }

  //single deme, constant N
  template< typename gamete_type,
	    typename glist_vector_type_allocator,
	    typename gamete_list_type_allocator,
	    typename mutation_list_type_allocator,
	    typename diploid_vector_type_allocator,
	    typename locus_vector_type_allocator,
	    typename diploid_fitness_function,
	    typename mutation_removal_policy,
	    typename mutation_model_container,
	    typename recombination_policy_container,
	    typename mutation_insertion_policy,
	    typename gamete_insertion_policy,
	    typename bw_locus_rec_fxn,
	    template<typename,typename> class gamete_list_type,
	    template<typename,typename> class glist_vector_type,
	    template<typename,typename> class mutation_list_type,
	    template<typename,typename> class diploid_vector_type,
	    template<typename,typename> class locus_vector_type>
  double
  sample_diploid(gsl_rng * r,
		 glist_vector_type< gamete_list_type<gamete_type,
		 gamete_list_type_allocator> ,
		 glist_vector_type_allocator > * gametes,
		 diploid_vector_type<locus_vector_type<std::pair<typename gamete_list_type< gamete_type,gamete_list_type_allocator >::iterator,
		 typename gamete_list_type< gamete_type,gamete_list_type_allocator >::iterator>,
		 locus_vector_type_allocator>,
		 diploid_vector_type_allocator> * diploids,
		 mutation_list_type<typename gamete_type::mutation_type,mutation_list_type_allocator > * mutations, 
		 const unsigned & N,
		 const double * mu,
		 const mutation_model_container & mmodel,
		 const recombination_policy_container & rec_policies,
		 const double * r_between_loci,
		 const bw_locus_rec_fxn & blrf,
		 const mutation_insertion_policy & mpolicy,
		 const gamete_insertion_policy & gpolicy_mut,
		 const diploid_fitness_function & ff,
		 const mutation_removal_policy & mp,
		 const double & f)
  {
    return sample_diploid(r,gametes,diploids,mutations,N,N,mu,mmodel,rec_policies,r_between_loci,blrf,mpolicy,gpolicy_mut,ff,mp,f);
  }
} //ns KTfwd
#endif
