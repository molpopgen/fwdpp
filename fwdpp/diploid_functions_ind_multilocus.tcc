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
	       const recombination_policy_container & rec_policies,
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
  
  typedef locus_vector_type<std::pair<typename gamete_list_type< gamete_type,gamete_list_type_allocator >::iterator,
				      typename gamete_list_type< gamete_type,gamete_list_type_allocator >::iterator>,
			    locus_vector_type_allocator> loci_ctr;

  typedef typename loci_ctr::iterator locus_itr;
  
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
 for( unsigned curr_dip = 0 ; curr_dip < N_next ; ++curr_dip )
   {
     assert(dptr==diploids->begin());
     assert( (dptr+curr_dip) < diploids->end() );
     size_t p1 = gsl_ran_discrete(r,lookup);
     size_t p2 = (gsl_rng_uniform(r) <= f) ? p1 : gsl_ran_discrete(r,lookup);
     assert(p1<parents.size());
     assert(p2<parents.size());
     
     //std::vector<unsigned> nrecs_p1( gametes->size() , 0u ),nrecs_p2( gametes->size,0u ); //store the number of recs per locus

     //Need to store a vector of the equivalent of p1g1,p1g2 out to png1,png2
     loci_ctr p1c( *(pptr+p1) ),
       p2c( *(pptr+p2) );

     /*
       Recombine -- updating via a bit field of 3 values.  
       The fields are:
       1.  Was last gamete g1 or g2?
       2.  Was last # crossovers even or odd?
       3.  Do we crossover b/w locus i and i-1?
     */
     int LW1=0,LW2=0;
     unsigned NR1=0,NR2=0;
     //Mendel:
     bool p1g1 = (gsl_rng_uniform(r) <= 0.5) ? true : false,
       p2g1 = (gsl_rng_uniform(r) <= 0.5) ? true : false;
     for( unsigned i = 0 ; i < p1c.size() ; ++i )
       {
	 locus_itr ptr2cdip = (dptr+curr_dip)->begin();
	 unsigned temp = rec_policies[i]( p1c[i].first, p1c[i].second );
	 //HARD PART IS HERE
	 if ( i > 1 )
	   {
	     int val = ( LW1 |= ( (NR1%2==0.) ? 2 : 0) ) |= ( (gsl_rng_uniform(r) <= *(r_between_loci+i-1)) ? 4 : 0 );
	     //int need_swap = (val != 2 && val != 4 && val != 7);
	     LW1 = (val != 2 && val != 4 && val != 7);
	     // if( LW1 )
	     //   {
	     // 	 //std::swap( p1c[i].first,p1c[i].second );
	     // 	 (ptr2cdip+i)->first = (p1g1) ? p1c[i].second : p1c[i].first;
	     //   }
	     // else
	     //   {
	     // 	 (ptr2cdip+i)->first = (p1g1) ? p1c[i].first : p1c[i].second;
	     //   }
	     
	     (ptr2cdip+i)->first = (LW1) ? ( (p1g1) ? p1c[i].second : p1c[i].first ) : ((p1g1) ? p1c[i].first : p1c[i].second);
	     //LW1 = need_swap;
	   }
	 NR1 = temp;
	 //(ptr2cdip+i)->first = (p1g1) ? p1c[i].first : p1c[i].second;

	 temp= rec_policies[i]( p2c[i].first, p2c[i].second );
	 //HARD PART IS HERE
	 if ( i > 1 )
	   {
	     int val = ( LW2 |= ( (NR2%2==0.) ? 2 : 0) ) |= ( (gsl_rng_uniform(r) <= *(r_between_loci+i-1)) ? 4 : 0 );
	     LW2 = (val != 2 && val != 4 && val != 7);
	     // if( LW2 )
	     //   {
	     // 	 //std::swap( p2c[i].first,p2c[i].second );
	     // 	 (ptr2cdip+i)->second = (p1g2) ? p2c[i].second : p2c[i].first;
	     //   }
	     // else
	     //   {
	     // 	 (ptr2cdip+i)->second = (p1g2) ? p2c[i].first : p2c[i].second;
	     //   }
	     (ptr2cdip+i)->second = (LW2) ? ((p1g2) ? p2c[i].second : p2c[i].first) : ((p1g2) ? p2c[i].first : p2c[i].second);
	     //LW2 = need_swap;
	   }
	 NR2 = temp;
	 //(ptr2cdip+i)->second = (p2g1) ? p2c[i].first : p2c[i].second;
	 
	 (ptr2cdip+i)->first->n++;
	 (ptr2cdip+i)->second->n++;
	 adjust_mutation_counts( (ptr2cdip+i)->first,1 );
	 adjust_mutation_counts( (ptr2cdip+i)->second,1 );
	 (ptr2cdip+i)->first = mutate_gamete( r,mu,(gametes+i),mutations,(ptr2cdip+i)->first,mmodel[i],mpolicy,gpolicy_mut);
	 (ptr2cdip+i)->second = mutate_gamete( r,mu,(gametes+i),mutations,(ptr2cdip+i)->second,mmodel[i],mpolicy,gpolicy_mut);
       }

     //single-locus version commented out beow  to use as reference
     // p1g1 = (pptr+p1)->first;
     // p1g2 = (pptr+p1)->second;
     // p2g1 = (pptr+p2)->first;
     // p2g2 = (pptr+p2)->second;
	
	// //Need a multilocus version of this bad boy -- or do we?
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
