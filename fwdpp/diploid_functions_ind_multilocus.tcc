//  -*- C++ -*- 
#ifndef __DIPLOID_FUNCTIONS_IND_MULTILOC_TCC__
#define __DIPLOID_FUNCTIONS_IND_MULTILOC_TCC__

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
	       const double & mu,
	       const mutation_model_container & mmodel,
	       const recombination_policy_container & rec_pol,
	       const double * r_between_loci,
	       const mutation_insertion_policy & mpolicy,
	       const gamete_insertion_policy & gpolicy_mut,
	       const diploid_fitness_function & ff,
	       const mutation_removal_policy & mp,
	       const double & f)
{
  assert(N_curr == diploids->size());
  
  for( typename mutation_list_type<typename gamete_type::mutation_type,mutation_list_type_allocator>::iterator itr = mutations->begin() ; 
       itr != mutations->end() ; ++itr )
    {
      itr->n = 0;
    }
  
  typedef glist_vector_type< gamete_list_type<gamete_type,gamete_list_type_allocator >,
			     glist_vector_type_allocator > multiloc_gcont;

typedef gamete_list_type<gamete_type,gamete_list_type_allocator > gamete_cont;

typedef diploid_vector_type<locus_vector_type<std::pair<typename gamete_list_type< gamete_type,gamete_list_type_allocator >::iterator,
							typename gamete_list_type< gamete_type,gamete_list_type_allocator >::iterator>,
					      locus_vector_type_allocator>,
			    diploid_vector_type_allocator> dipctr;
 
 typedef typename locus_vector_type<std::pair<typename gamete_list_type< gamete_type,gamete_list_type_allocator >::iterator,
					      typename gamete_list_type< gamete_type,gamete_list_type_allocator >::iterator>,
				    locus_vector_type_allocator>::iterator locus_itr;
 
 std::vector<double> fitnesses(diploids->size());
 double wbar = 0.;
 
 typename dipctr::iterator dptr = diploids->begin();
 for( unsigned i = 0 ; i < N_curr ; ++i,++dptr )
   {
     for( locus_itr j = dptr->begin() ; j != dptr->end() ; ++j )
       {
	 j->first->n = 0;
	 j->second->n = 0;
       }
     fitnesses[i] = ff( dptr );
     wbar += fitnesses[i];
   }
 wbar /= double(diploids->size());
 dptr = diploids->begin();
 
#ifndef NDEBUG
 for ( typename multiloc_gcont::iterator i = gametes->begin() ; i != gametes->end() ; ++i )
   {
     for(typename gamete_list_type<gamete_type,gamete_list_type_allocator >::iterator itr = i->begin() ; itr != i->end() ;++itr)
       {
	 assert(itr->n==0);
       }
   }
#endif
 
 gsl_ran_discrete_t * lookup = gsl_ran_discrete_preproc(fitnesses.size(),&fitnesses[0]);
 
 dipctr parents(*diploids); //copy the parents
 const typename dipctr::iterator pptr = parents.begin();
 
 //Change the population size
 if( diploids->size() != N_next)
   {
     diploids->resize(N_next);
     dptr = diploids->begin();
   }
 unsigned NREC=0;
 assert(diploids->size()==N_next);
 typename gamete_list_type< gamete_type,gamete_list_type_allocator >::iterator p1g1,p1g2,p2g1,p2g2;
 for( unsigned i = 0 ; i < N_next ; ++i )
   {
     assert(dptr==diploids->begin());
     assert( (dptr+i) < diploids->end() );
     size_t p1 = gsl_ran_discrete(r,lookup);
     size_t p2 = (gsl_rng_uniform(r) <= f) ? p1 : gsl_ran_discrete(r,lookup);
     assert(p1<parents.size());
     assert(p2<parents.size());
     
     //single-locus version commented out beow
     // p1g1 = (pptr+p1)->first;
     // p1g2 = (pptr+p1)->second;
     // p2g1 = (pptr+p2)->first;
     // p2g2 = (pptr+p2)->second;
	
	// //Need a multilocus version of this bad boy
	// NREC += rec_pol(p1g1,p1g2);
	// NREC += rec_pol(p2g1,p2g2);
	
	// (dptr+i)->first = (gsl_rng_uniform(r) <= 0.5) ? p1g1 : p1g2;
	// (dptr+i)->second = (gsl_rng_uniform(r) <= 0.5) ? p2g1 : p2g2;
	
	// (dptr+i)->first->n++;
	// assert( (dptr+i)->first->n > 0 );
	// assert( (dptr+i)->first->n <= 2*N_next );
	// (dptr+i)->second->n++;
	// assert( (dptr+i)->second->n > 0 );
	// assert( (dptr+i)->second->n <= 2*N_next );
	
	// adjust_mutation_counts((dptr+i)->first,1);
	// adjust_mutation_counts((dptr+i)->second,1);
	
	// //now, add new mutations
	// (dptr+i)->first = mutate_gamete(r,mu,gametes,mutations,(dptr+i)->first,mmodel,mpolicy,gpolicy_mut);
	// (dptr+i)->second = mutate_gamete(r,mu,gametes,mutations,(dptr+i)->second,mmodel,mpolicy,gpolicy_mut);
      }
#ifndef NDEBUG
    for( unsigned i = 0 ; i < diploids->size() ; ++i )
      {
	//check that gamete counts are all ok
      }
#endif

    for ( typename multiloc_gcont::iterator i = gametes->begin() ; i != gametes->end() ; ++i )
      {
       	i->remove_if(boost::bind(n_is_zero(),_1));
	for (typename gamete_cont::iterator gptr = i->begin() ; gptr != i->end() ; ++gptr )
	  {
	    gptr->mutations.erase( std::remove_if(gptr->mutations.begin(),gptr->mutations.end(),mp), gptr->mutations.end() );
	    gptr->smutations.erase( std::remove_if(gptr->smutations.begin(),gptr->smutations.end(),mp), gptr->smutations.end() );
	  }
      }
    gsl_ran_discrete_free(lookup);
    return wbar;
  }

} //ns KTfwd
#endif
