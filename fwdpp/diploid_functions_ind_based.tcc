//  -*- C++ -*- 
#ifndef __DIPLOID_FUNCTIONS_IND_BASED_TCC__
#define __DIPLOID_FUNCTIONS_IND_BASED_TCC__

#include <fwdpp/internal/gsl_discrete.hpp>
#include <fwdpp/internal/diploid_fitness_dispatch.hpp>
#include <fwdpp/internal/gamete_lookup_table.hpp>
#include <fwdpp/internal/gamete_cleaner.hpp>

namespace KTfwd
{
  //single deme, constant N
  template< typename gamete_type,
	    typename gamete_list_type_allocator,
	    typename mutation_list_type_allocator,
	    typename diploid_geno_t,
	    typename diploid_vector_type_allocator,
	    typename diploid_fitness_function,
	    typename mutation_removal_policy,
	    typename mutation_model,
	    typename recombination_policy,
	    typename mutation_insertion_policy,
	    typename gamete_insertion_policy,
	    template<typename,typename> class gamete_list_type,
	    template<typename,typename> class mutation_list_type,
	    template<typename,typename> class diploid_vector_type>
  double
  sample_diploid(gsl_rng * r,
		 gamete_list_type<gamete_type,gamete_list_type_allocator > * gametes,
		 diploid_vector_type<diploid_geno_t,diploid_vector_type_allocator> * diploids,
		 mutation_list_type<typename gamete_type::mutation_type,mutation_list_type_allocator > * mutations, 
		 const unsigned & N_curr, 
		 const double & mu,
		 const mutation_model & mmodel,
		 const recombination_policy & rec_pol,
		 const mutation_insertion_policy & mpolicy,
		 const gamete_insertion_policy & gpolicy_mut,
		 const diploid_fitness_function & ff,
		 const mutation_removal_policy & mp,
		 const double & f)
  {
    //run changing N version with N_next == N_curr
    return sample_diploid(r,gametes,diploids,mutations,N_curr,N_curr,mu,mmodel,rec_pol,mpolicy,
			  gpolicy_mut,ff,mp,f);
  }

  //single deme, N changing
  template< typename gamete_type,
	    typename gamete_list_type_allocator,
	    typename mutation_list_type_allocator,
	    typename diploid_geno_t,
	    typename diploid_vector_type_allocator,
	    typename diploid_fitness_function,
	    typename mutation_removal_policy,
	    typename mutation_model,
	    typename recombination_policy,
	    typename mutation_insertion_policy,
	    typename gamete_insertion_policy,
	    template<typename,typename> class gamete_list_type,
	    template<typename,typename> class mutation_list_type,
	    template<typename,typename> class diploid_vector_type>
  double
  sample_diploid(gsl_rng * r,
		 gamete_list_type<gamete_type,gamete_list_type_allocator > * gametes,
		 diploid_vector_type<diploid_geno_t,diploid_vector_type_allocator> * diploids,
		 mutation_list_type<typename gamete_type::mutation_type,mutation_list_type_allocator > * mutations, 
		 const unsigned & N_curr, 
		 const unsigned & N_next, 
		 const double & mu,
		 const mutation_model & mmodel,
		 const recombination_policy & rec_pol,
		 const mutation_insertion_policy & mpolicy,
		 const gamete_insertion_policy & gpolicy_mut,
		 const diploid_fitness_function & ff,
		 const mutation_removal_policy & mp,
		 const double & f)
  {
    assert(N_curr == diploids->size());

    std::vector<double> fitnesses(diploids->size());
    double wbar = 0.;
    
    auto dptr = diploids->begin();
    for( unsigned i = 0 ; i < N_curr ; ++i )
      {
	(dptr+i)->first->n = 0;
	(dptr+i)->second->n = 0;
	fitnesses[i] = fwdpp_internal::diploid_fitness_dispatch(ff,(dptr+i),
								typename traits::is_custom_diploid_t<diploid_geno_t>::type());
	wbar += fitnesses[i];
      }
    wbar /= double(diploids->size());

#ifndef NDEBUG
    std::for_each(gametes->cbegin(),gametes->cend(),[](decltype((*gametes->cbegin())) __g) {
	assert( !__g.n ); } );
#endif
    fwdpp_internal::gsl_ran_discrete_t_ptr lookup(gsl_ran_discrete_preproc(fitnesses.size(),&fitnesses[0]));
    auto parents(*diploids); //copy the parents
    auto pptr = parents.begin();
    
    //Change the population size
    if( diploids->size() != N_next)
      {
	diploids->resize(N_next);
	dptr = diploids->begin();
      }
    unsigned NREC=0;
    assert(diploids->size()==N_next);
    decltype( gametes->begin() ) p1g1,p1g2,p2g1,p2g2;

    auto gamete_lookup = fwdpp_internal::gamete_lookup_table(gametes);

    for( unsigned i = 0 ; i < N_next ; ++i )
      {
	assert(dptr==diploids->begin());
	assert( (dptr+i) < diploids->end() );
	size_t p1 = gsl_ran_discrete(r,lookup.get());
#ifdef FWDPP_COMPAT_0_3_0
	size_t p2 = (gsl_rng_uniform(r) <= f) ? p1 : gsl_ran_discrete(r,lookup.get());
#else
	size_t p2 = (f==1. || (f>0. && gsl_rng_uniform(r)<=f)) ? p1 : gsl_ran_discrete(r,lookup.get());
#endif
	assert(p1<parents.size());
	assert(p2<parents.size());
	
	p1g1 = (pptr+typename decltype(pptr)::difference_type(p1))->first;
	p1g2 = (pptr+typename decltype(pptr)::difference_type(p1))->second;
	p2g1 = (pptr+typename decltype(pptr)::difference_type(p2))->first;
	p2g2 = (pptr+typename decltype(pptr)::difference_type(p2))->second;

	if(gsl_rng_uniform(r)<=0.5) std::swap(p1g1,p1g2);
	if(gsl_rng_uniform(r)<=0.5) std::swap(p2g1,p2g2);

	NREC += rec_pol(p1g1,p1g2,gamete_lookup);
	NREC += rec_pol(p2g1,p2g2,gamete_lookup);

	(dptr+i)->first = p1g1;
	(dptr+i)->second = p2g1;
	
	(dptr+i)->first->n++;
	assert( (dptr+i)->first->n > 0 );
	assert( (dptr+i)->first->n <= 2*N_next );
	(dptr+i)->second->n++;
	assert( (dptr+i)->second->n > 0 );
	assert( (dptr+i)->second->n <= 2*N_next );

	//now, add new mutations
	(dptr+i)->first = mutate_gamete(r,mu,gametes,mutations,(dptr+i)->first,mmodel,mpolicy,gpolicy_mut);
	(dptr+i)->second = mutate_gamete(r,mu,gametes,mutations,(dptr+i)->second,mmodel,mpolicy,gpolicy_mut);
      }
#ifndef NDEBUG
    for( unsigned i = 0 ; i < diploids->size() ; ++i )
      {
	assert( (dptr+i)->first->n > 0 );
	assert( (dptr+i)->first->n <= 2*N_next );
	assert( (dptr+i)->second->n > 0 );
	assert( (dptr+i)->second->n <= 2*N_next );
      }
    for( auto itr = mutations->begin() ; itr != mutations->end() ; ++itr ) assert( itr->n <= 2*N_next );
#endif
    decltype(gametes->begin()) temp;

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
    assert(check_sum(gametes,2*N_next));
    return wbar;
  }
  
  //Metapopulation version of sample_diploid for individual-based simulations and constant N
  template< typename gamete_type,
	    typename metapop_diploid_vector_type_allocator,
	    typename gamete_list_type_allocator,
	    typename mutation_list_type_allocator,
	    typename diploid_geno_t,
	    typename diploid_vector_type_allocator,
	    typename diploid_fitness_function_container,
	    typename mutation_removal_policy,
	    typename mutation_model,
	    typename recombination_policy,
	    typename migration_policy,
	    typename mutation_insertion_policy,
	    typename gamete_insertion_policy,
	    template<typename,typename> class gamete_list_type,
	    template<typename,typename> class mutation_list_type,
	    template<typename,typename> class diploid_vector_type,
	    template<typename,typename> class metapop_diploid_vector_type>
  std::vector< double >
  sample_diploid(gsl_rng * r,
		 gamete_list_type<gamete_type,gamete_list_type_allocator> * metapop,
		 metapop_diploid_vector_type < diploid_vector_type<diploid_geno_t,diploid_vector_type_allocator>,metapop_diploid_vector_type_allocator > * diploids,
		 mutation_list_type<typename gamete_type::mutation_type,mutation_list_type_allocator > * mutations, 
		 const unsigned * N_curr, 
		 const double & mu,
		 const mutation_model & mmodel,
		 const recombination_policy & rec_pol,
		 const mutation_insertion_policy & mpolicy,
		 const gamete_insertion_policy & gpolicy_mut,
		 const diploid_fitness_function_container & ffs,
		 const mutation_removal_policy & mp,
		 const migration_policy & mig,
		 const double * f)
  {
    //run changing-N version with no change in N
    return sample_diploid(r,metapop,diploids,mutations,N_curr,N_curr,mu,mmodel,rec_pol,mpolicy,
			  gpolicy_mut,ffs,mp,mig,f);
  }
  
  //Metapopulation version of sample_diploid for individual-based simulations with changing population size
  template< typename gamete_type,
	    typename metapop_diploid_vector_type_allocator,
	    typename gamete_list_type_allocator,
	    typename mutation_list_type_allocator,
	    typename diploid_geno_t,
	    typename diploid_vector_type_allocator,
	    typename diploid_fitness_function_container,
	    typename mutation_removal_policy,
	    typename mutation_model,
	    typename recombination_policy,
	    typename migration_policy,
	    typename mutation_insertion_policy,
	    typename gamete_insertion_policy,
	    template<typename,typename> class gamete_list_type,
	    template<typename,typename> class mutation_list_type,
	    template<typename,typename> class diploid_vector_type,
	    template<typename,typename> class metapop_diploid_vector_type>
  std::vector< double >
  sample_diploid(gsl_rng * r,
		 gamete_list_type<gamete_type,gamete_list_type_allocator> * metapop,
		 metapop_diploid_vector_type < diploid_vector_type<diploid_geno_t, diploid_vector_type_allocator>,metapop_diploid_vector_type_allocator > * diploids,
		 mutation_list_type<typename gamete_type::mutation_type,mutation_list_type_allocator > * mutations, 
		 const unsigned * N_curr, 
		 const unsigned * N_next, 
		 const double & mu,
		 const mutation_model & mmodel,
		 const recombination_policy & rec_pol,
		 const mutation_insertion_policy & mpolicy,
		 const gamete_insertion_policy & gpolicy_mut,
		 const diploid_fitness_function_container & ffs,
		 const mutation_removal_policy & mp,
		 const migration_policy & mig,
		 const double * f)
	    {
	      //get the fitnesses for each diploid in each deme and make the lookup table of parental fitnesses
	      using lookup_t = fwdpp_internal::gsl_ran_discrete_t_ptr;
	      std::vector<lookup_t> lookups;
	      std::vector<double> wbars(diploids->size(),0);
	      typename decltype(diploids->begin())::difference_type popindex = 0;

	      //get max N
	      unsigned mN=0;
	      for( unsigned i=0;i<diploids->size();++i )
		{
		  if( *(N_curr+i) > mN )
		    {
		      mN = *(N_curr+i);
		    }
		}
	      double * fitnesses = new double[mN];
	      
	      for( auto dptr = diploids->begin() ; dptr != diploids->end() ; ++dptr, ++popindex )
		{
		  unsigned demesize = *(N_curr+popindex);
		  assert( demesize == dptr->size() );
		  size_t ith_dip = 0;
		  for( auto gptr = dptr->begin() ; 
		       gptr != dptr->end() ; ++gptr,++ith_dip )
		    {
		      fitnesses[ith_dip] = fwdpp_internal::diploid_fitness_dispatch(ffs[typename diploid_fitness_function_container::size_type(popindex)],gptr,
										    typename traits::is_custom_diploid_t<diploid_geno_t>::type());
		      wbars[std::vector<double>::size_type(popindex)]+=fitnesses[ith_dip];
		      gptr->first->n = 0;
		      gptr->second->n = 0;
		    }
		  wbars[std::vector<double>::size_type(popindex)] /= double( demesize );
		  lookups.emplace_back( lookup_t(gsl_ran_discrete_preproc(demesize,fitnesses)) );
		}
	      delete [] fitnesses;
	      
	      assert(lookups.size() == diploids->size());
	      //copy diploids into temporary parents
	      auto parents(*diploids);
	      
	      //Update the diploids
	      popindex = 0;
	      unsigned NREC=0;

	      decltype(metapop->begin()) p1g1,p1g2,p2g1,p2g2;
	      auto gamete_lookup = fwdpp_internal::gamete_lookup_table(metapop);
	      for( auto ptr = diploids->begin() ; ptr != diploids->end() ; ++ptr,++popindex )
		{
		  unsigned demesize =*(N_next+popindex);
		  if( demesize != *(N_curr+popindex) )
		    {
		      ptr->resize(demesize);
		    }
		  auto dptr = ptr->begin();

		  for( unsigned i = 0 ; i < demesize ; ++i )
		    {
		      /* Figure out if parent 1 is migrant or not.

			A migration policy takes the current deme (popindex) as
			an argument.  It returns popindex if there is no migration,
			else it returns the index of the deme of a migrant parent
		      */
		      decltype(popindex) deme_first_parent = decltype(popindex)(mig(size_t(popindex))),deme_other_parent=popindex;
		      auto pptr=(parents.begin()+typename decltype(parents.begin())::difference_type(deme_first_parent))->begin();
		      typename decltype(pptr)::difference_type p1 = 
			typename decltype(pptr)::difference_type(gsl_ran_discrete(r,lookups[std::vector<lookup_t>::size_type(deme_first_parent)].get())),p2;

		      p1g1 = (pptr+p1)->first;
		      p1g2 = (pptr+p1)->second;

		      /*
			If the individual is not inbred, then we pick a 
			deme from the migration policy for parent 2
		      */
		      auto pptr2=(parents.begin()+typename decltype(parents.begin())::difference_type(deme_other_parent))->end();
#ifdef FWDPP_COMPAT_0_3_0
		      if( f != nullptr && gsl_rng_uniform(r) <= *(f + popindex ) ) //individual is inbred
#else
			if( f != nullptr && ( *(f + popindex)==1. || (*(f + popindex)>0. && gsl_rng_uniform(r) <= *(f + popindex)) ) ) //individual is inbred
#endif
			  {
			    pptr2=(parents.begin()+typename decltype(parents.begin())::difference_type(popindex))->begin();
			    p2=p1;
			  }
		      else
			{
			  deme_other_parent = decltype(deme_other_parent)(mig(size_t(popindex)));
			  assert(deme_other_parent>=0);
			  assert( decltype(diploids->size())(deme_other_parent) < diploids->size() );
			  pptr2 = (parents.begin() + deme_other_parent)->begin();
			  p2 = decltype(p2)(gsl_ran_discrete(r,lookups[std::vector<lookup_t>::size_type(deme_other_parent)].get()));
			  assert( (pptr2+p2) < (parents.begin() + typename decltype(parents.begin())::difference_type(deme_other_parent))->end() );
			}
		      assert( pptr2 != (parents.begin() + typename decltype(parents.begin())::difference_type(deme_other_parent))->end() );

		      p2g1 = (pptr2+p2)->first;
		      p2g2 = (pptr2+p2)->second;

		      //0.3.3: Do "Mendel" now...
		      if(gsl_rng_uniform(r)<=0.5) std::swap(p1g1,p1g2);
		      if(gsl_rng_uniform(r)<=0.5) std::swap(p2g1,p2g2);
		      
		      NREC += rec_pol(p1g1,p1g2,gamete_lookup);
		      NREC += rec_pol(p2g1,p2g2,gamete_lookup);

		      (dptr+i)->first = p1g1;
		      (dptr+i)->second = p2g1;
		      assert( std::find( (metapop)->begin(), (metapop)->end(), *( (dptr+i)->second ) )
			      != (metapop)->end() );

		      (dptr+i)->first->n++;
		      (dptr+i)->second->n++;

		      (dptr+i)->first = mutate_gamete(r,mu,metapop,mutations,(dptr+i)->first,mmodel,mpolicy,gpolicy_mut);
		      (dptr+i)->second = mutate_gamete(r,mu,metapop,mutations,(dptr+i)->second,mmodel,mpolicy,gpolicy_mut);
		    }
		}
	      for( auto itr = metapop->begin() ; itr != metapop->end() ; )
		{
		  if(!itr->n) itr = metapop->erase(itr);
		  else
		    {
		      adjust_mutation_counts(itr,itr->n);
		      ++itr; 
		    }
		}
	      fwdpp_internal::gamete_cleaner(metapop,mp,typename std::is_same<mutation_removal_policy,KTfwd::remove_nothing >::type());
	      return wbars;
	    }
}

#endif
  
