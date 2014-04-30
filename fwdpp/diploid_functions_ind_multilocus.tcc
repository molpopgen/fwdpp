//  -*- C++ -*- 
#ifndef __DIPLOID_FUNCTIONS_IND_MULTILOC_TCC__
#define __DIPLOID_FUNCTIONS_IND_MULTILOC_TCC__

#include <fwdpp/mutation.hpp>

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
	       const double * mu,
	       const mutation_model_container & mmodel,
	       const recombination_policy_container & rec_policies,
	       const double * r_between_loci,
	       const mutation_insertion_policy & mpolicy,
	       const gamete_insertion_policy & gpolicy_mut,
	       const diploid_fitness_function & ff,
	       const mutation_removal_policy & mp,
	       const double & f)
{
  typedef gamete_list_type<gamete_type,gamete_list_type_allocator > gamete_cont;
  
  typedef glist_vector_type< gamete_cont, glist_vector_type_allocator > multiloc_gcont;

  typedef typename gamete_list_type< gamete_type,gamete_list_type_allocator >::iterator glist_itr;

  typedef locus_vector_type<std::pair<glist_itr,glist_itr>,
			    locus_vector_type_allocator> loci_ctr;

  typedef diploid_vector_type<loci_ctr,
			      diploid_vector_type_allocator> dipctr;
 
  typedef typename loci_ctr::iterator locus_itr;

#ifndef NDEBUG
  assert(diploids->size()==N_curr);
#endif

  //Adjust all mutation frequencies down to 0
  for( typename mutation_list_type<typename gamete_type::mutation_type,mutation_list_type_allocator>::iterator itr = mutations->begin() ; 
       itr != mutations->end() ; ++itr )
    {
      itr->n = 0;
    }
  
  //Vector of parental fitnesses
  std::vector<double> fitnesses(N_curr);
  double wbar = 0.;
  
  //Go over parents
  typename dipctr::iterator dptr = diploids->begin();
  for( unsigned i = 0 ; i < N_curr ; ++i,++dptr )
  {
    //set parental gamete counts to 0 for each locus
    for( locus_itr j = dptr->begin() ; j != dptr->end() ; ++j )
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
  for ( typename multiloc_gcont::const_iterator i = gametes->begin() ; i != gametes->end() ; ++i )
    {
      for (typename gamete_cont::const_iterator gptr = i->begin() ; gptr != i->end() ; ++gptr )
	{
	  assert( ! gptr->n );
	}
    }
#endif

  wbar /= double(diploids->size());
  dptr = diploids->begin();  //reset to beginning of diploids
  
  gsl_ran_discrete_t * lookup = gsl_ran_discrete_preproc(fitnesses.size(),&fitnesses[0]);
 
  dipctr parents(*diploids); //Copy the parents.  Exact copy of diploids--same fitnesses, etc.
  const typename dipctr::iterator pptr = parents.begin(); //the first parent
 
  //Change the population size and reset dptr to avoid iterator invalidation
  if( diploids->size() != N_next)
    {
      diploids->resize(N_next);
      dptr = diploids->begin();
    }
  
  assert(diploids->size()==N_next);
 
  //WTF is this for?
  /*
  std::vector< gamete_cont * > ptr_to_gametes;
  for( typename multiloc_gcont::iterator itr = gametes->begin() ; itr != gametes->end() ; ++itr )
    {
      ptr_to_gametes.push_back( &*itr );
   }
  */
 for( unsigned curr_dip = 0 ; curr_dip < N_next ; ++curr_dip )
   {
     assert(dptr==diploids->begin());
     assert( (dptr+curr_dip) < diploids->end() );
     size_t p1 = gsl_ran_discrete(r,lookup);
     size_t p2 = (gsl_rng_uniform(r) <= f) ? p1 : gsl_ran_discrete(r,lookup);
     assert(p1<N_curr);
     assert(p2<N_curr);
     
     //Need to store a vector of the equivalent of p1g1,p1g2 out to p1gn,p2gn
     //This is a trivial copying of iterators, so not that expensive
     loci_ctr p1c( *(pptr+p1) ),
       p2c( *(pptr+p2) );
     assert( p1c == *(pptr+p1) );
     assert( p2c == *(pptr+p2) );

     //Using just the routines below give correct E[S] for n = 20.  It must be the bitset nonsense that is incorrect.
     //Temporary debug mode: drift and mutation only!!!!  OK--this seems good based on some limited testing
     //recombination, too! OK--this seems good based on some limited testing.
     //All testing so far based on n=2 using strobeck_morgan.cc

     //For 3 locus (test/strobeck_morgan.cc), the dist. of 
     //summary stats looks great vis-a-vis ms for N=20,theta=20,rho=1, N=1e3,10N gens
     bool p1g1 = (gsl_rng_uniform(r) <= 0.5) ? true : false,
       p2g1 = (gsl_rng_uniform(r) <= 0.5) ? true : false;
     locus_itr ptr2cdip = (dptr+curr_dip)->begin();
     bool LO1 = true, LO2 = true;
     for ( unsigned i = 0 ; i < p1c.size() ; ++i )
       {
	 unsigned temp = rec_policies[i]( p1c[i].first, p1c[i].second );
	 if ( i > 0 )
	   {
	     unsigned nrbw = gsl_ran_poisson(r,r_between_loci[i-1]);
	     bool obw = (nrbw%2!=0) ? true : false;
	     //bool dummy = p1g1;
	     p1g1 = (LO1) ? !p1g1 : p1g1;
	     p1g1 = (obw) ? !p1g1 : p1g1;
	     /*
	     if ( nrbw )
	       {
		 std::cerr << nrbw << ' ' << obw << ' ' << dummy << ' ' << p1g1 << '\n';
	       }
	     */
	   }
	     /*(ptr2cdip+i)->first = (p1g1) ? p1c[i].first : p1c[i].second;
	   }
	 else
	   {
	     (ptr2cdip+i)->first = (p1g1) ? p1c[i].first : p1c[i].second;
	     }*/
	 (ptr2cdip+i)->first = (p1g1) ? p1c[i].first : p1c[i].second;
	 LO1 = (temp % 2 != 0.) ? true : false;

	 temp = rec_policies[i]( p2c[i].first, p2c[i].second );
	 if ( i > 0 )
	   {
	     unsigned nrbw = gsl_ran_poisson(r,r_between_loci[i-1]);
	     bool obw = (nrbw%2!=0) ? true : false;
	     p2g1 = (LO1) ? !p2g1 : p2g1;
	     p2g1 = (obw) ? !p2g1 : p2g1;
	   }
	 (ptr2cdip+i)->second = (p2g1) ? p2c[i].first : p2c[i].second;

	 (ptr2cdip+i)->first->n++;
	 (ptr2cdip+i)->second->n++;

	 adjust_mutation_counts( (ptr2cdip+i)->first,1 );
	 adjust_mutation_counts( (ptr2cdip+i)->second,1 );

	 (ptr2cdip+i)->first = mutate_gamete( r,mu[i],&*(gametes->begin()+i),mutations,(ptr2cdip+i)->first,mmodel[i],mpolicy,gpolicy_mut);
	 (ptr2cdip+i)->second = mutate_gamete( r,mu[i],&*(gametes->begin()+i),mutations,(ptr2cdip+i)->second,mmodel[i],mpolicy,gpolicy_mut);	 
       }
     /*
       Recombine -- updating via a bit field of 3 values.  
       The fields are:
       1.  Was last gamete pulled from p1c or p2c?  This is LW1/2, which takes on value 0 or 1 (false/true)
       2.  Was last # crossovers even or odd?  This is NR1/2, which takes on value 0 or 2 (false/true)
       3.  Do we crossover b/w locus i and i-1?   This is determined by sampling from r_between_loci and takes on value 0 or 4 (false/true)

       The above establises a bitset that may take on values 0 thru 7.
     */
     /*
     int LW1=0,LW2=0;
     unsigned NR1=0,NR2=0;
     //Mendel:
     bool p1g1 = (gsl_rng_uniform(r) <= 0.5) ? true : false,
       p2g1 = (gsl_rng_uniform(r) <= 0.5) ? true : false;
     locus_itr ptr2cdip = (dptr+curr_dip)->begin();
     for( unsigned i = 0 ; i < p1c.size() ; ++i )
       {
	 unsigned temp = rec_policies[i]( p1c[i].first, p1c[i].second );
	 if ( i > 0 )
	   {
	     int val = ( LW1 |= ( (NR1%2==0.) ? 2 : 0) ) |= ( (gsl_rng_uniform(r) <= *(r_between_loci+i-1)) ? 4 : 0 );
	     LW1 = (val != 2 && val != 4 && val != 7) ? 1 : 0;
	     (ptr2cdip+i)->first = (LW1) ? ( (p1g1) ? p1c[i].second : p1c[i].first ) : ((p1g1) ? p1c[i].first : p1c[i].second);
	   }
	 else
	   {
	     (ptr2cdip+i)->first = (p1g1) ? p1c[i].first : p1c[i].second;
	     LW1 = int(!p1g1);
	   }
	 NR1 = temp;

	 temp = rec_policies[i]( p2c[i].first, p2c[i].second );
	 if ( i > 0 )
	   {
	     int val = ( LW2 |= ( (NR2%2==0.) ? 2 : 0) ) |= ( (gsl_rng_uniform(r) <= *(r_between_loci+i-1)) ? 4 : 0 );
	     LW2 = (val != 2 && val != 4 && val != 7) ? 1 : 0;
	     (ptr2cdip+i)->second = (LW2) ? ((p2g1) ? p2c[i].second : p2c[i].first) : ((p2g1) ? p2c[i].first : p2c[i].second);
	   }
	 else
	   {
	     (ptr2cdip+i)->second = ( p2g1 ) ? p2c[i].first : p2c[i].second;
	     LW2 = int(!p2g1);
	   }
	 NR2 = temp;
	 
	 (ptr2cdip+i)->first->n++;
	 (ptr2cdip+i)->second->n++;

	 adjust_mutation_counts( (ptr2cdip+i)->first,1 );
	 adjust_mutation_counts( (ptr2cdip+i)->second,1 );

	 //(ptr2cdip+i)->first = mutate_gamete( r,mu[i],ptr_to_gametes[i],mutations,(ptr2cdip+i)->first,mmodel[i],mpolicy,gpolicy_mut);
	 //(ptr2cdip+i)->second = mutate_gamete( r,mu[i],ptr_to_gametes[i],mutations,(ptr2cdip+i)->second,mmodel[i],mpolicy,gpolicy_mut);
	 (ptr2cdip+i)->first = mutate_gamete( r,mu[i],&*(gametes->begin()+i),mutations,(ptr2cdip+i)->first,mmodel[i],mpolicy,gpolicy_mut);
	 (ptr2cdip+i)->second = mutate_gamete( r,mu[i],&*(gametes->begin()+i),mutations,(ptr2cdip+i)->second,mmodel[i],mpolicy,gpolicy_mut);
       }
     */
   }
    for ( typename multiloc_gcont::iterator i = gametes->begin() ; i != gametes->end() ; ++i )
      {
       	i->remove_if(boost::bind(n_is_zero(),_1));
	for (typename gamete_cont::iterator gptr = i->begin() ; gptr != i->end() ; ++gptr )
	  {
	    gptr->mutations.erase( std::remove_if(gptr->mutations.begin(),gptr->mutations.end(),mp), gptr->mutations.end() );
	    gptr->smutations.erase( std::remove_if(gptr->smutations.begin(),gptr->smutations.end(),mp), gptr->smutations.end() );
	  }
      }
#ifndef NDEBUG
    for ( typename multiloc_gcont::const_iterator i = gametes->begin() ; i != gametes->end() ; ++i )
      {
	unsigned sum = 0;
	for (typename gamete_cont::const_iterator gptr = i->begin() ; gptr != i->end() ; ++gptr )
	  {
	    //make sure that mutation frequencies are >= gamete frequencies
	    for( unsigned j = 0 ; j < gptr->mutations.size() ; ++j )
	      {
		assert( gptr->mutations[j]->n >= gptr->n );
	      }
	    for( unsigned j = 0 ; j < gptr->smutations.size() ; ++j )
	      {
		assert( gptr->smutations[j]->n >= gptr->n );
	      }
	    
	    sum += gptr->n;
	    assert( sum && sum <= 2*N_next );
	  }
	assert(sum == 2*N_next);
      }
#endif
    gsl_ran_discrete_free(lookup);
    return wbar;
  }

} //ns KTfwd
#endif
