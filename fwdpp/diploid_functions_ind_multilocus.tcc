//  -*- C++ -*- 
#ifndef __DIPLOID_FUNCTIONS_IND_MULTILOC_TCC__
#define __DIPLOID_FUNCTIONS_IND_MULTILOC_TCC__

#include <fwdpp/mutation.hpp>
#include <fwdpp/internal/gsl_discrete.hpp>
#include <fwdpp/internal/multilocus_rec.hpp>
#include <fwdpp/internal/gamete_lookup_table.hpp>
#include <fwdpp/internal/gamete_cleaner.hpp>

namespace KTfwd
{
  //single deme, N changing
  template< typename diploid_geno_t,
	    typename gamete_type,
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
	    template<typename,typename> class mutation_list_type,
	    template<typename,typename> class diploid_vector_type,
	    template<typename,typename> class locus_vector_type>
  double
  sample_diploid(gsl_rng * r,
		 gamete_list_type<gamete_type, gamete_list_type_allocator> * gametes,
		 diploid_vector_type<locus_vector_type<diploid_geno_t,locus_vector_type_allocator>,diploid_vector_type_allocator> * diploids,
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
	fitnesses[i] += ff( dptr );
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
	    assert( ! i->n );
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

    auto gamete_lookup = fwdpp_internal::gamete_lookup_table(gametes);
    for( unsigned curr_dip = 0 ; curr_dip < N_next ; ++curr_dip )
      {
	assert(dptr==diploids->begin());
	assert( (dptr+curr_dip) < diploids->end() );

	//Choose the two parents
	typename decltype(pptr)::difference_type p1 = decltype(p1)(gsl_ran_discrete(r,lookup.get()));
#ifdef FWDPP_COMPAT_0_3_0
	decltype(p1) p2  = (gsl_rng_uniform(r) <= f) ? p1 : decltype(p1)(gsl_ran_discrete(r,lookup.get()));
#else
	decltype(p1) p2 = (f==1. || (f>0. && gsl_rng_uniform(r)<=f)) ? p1 :  decltype(p1)(gsl_ran_discrete(r,lookup.get()));
#endif
	assert(p1<N_curr);
	assert(p2<N_curr);
     
	auto cdip = (dptr+curr_dip);
	fwdpp_internal::multilocus_rec_mut(r,*(pptr+p1),*(pptr+p2),cdip,gamete_lookup,
					   rec_policies,blrf,r_between_loci,
					   ((gsl_rng_uniform(r)<=0.5)?1:0),
					   ((gsl_rng_uniform(r)<=0.5)?1:0),
					   gametes,mutations,mu,mmodel,mpolicy,gpolicy_mut
					   );
      }

    //0.3.3: simpler!
    for( auto itr = gametes->begin() ; itr != gametes->end() ; )
      {
	if(!itr->n) itr = gametes->erase(itr);
	else
	  {
	    adjust_mutation_counts(itr,itr->n);
	    ++itr; 
	  }
      }
    fwdpp_internal::gamete_cleaner(gametes,mp,typename std::is_same<mutation_removal_policy,KTfwd::remove_nothing >::type());
    return wbar;
  }

  //single deme, constant N
  template< typename diploid_geno_t,
	    typename gamete_type,
	    //typename glist_vector_type_allocator,
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
	    //template<typename,typename> class glist_vector_type,
	    template<typename,typename> class mutation_list_type,
	    template<typename,typename> class diploid_vector_type,
	    template<typename,typename> class locus_vector_type>
  double
  sample_diploid(gsl_rng * r,
		 //IDEA:
		 //glist_vector_type< gamete_list_type<gamete_type,
		 //gamete_list_type_allocator> ,
		 //glist_vector_type_allocator > * gametes,
		 gamete_list_type<gamete_type, gamete_list_type_allocator> * gametes,
		 diploid_vector_type<locus_vector_type<diploid_geno_t,locus_vector_type_allocator>,diploid_vector_type_allocator> * diploids,
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
