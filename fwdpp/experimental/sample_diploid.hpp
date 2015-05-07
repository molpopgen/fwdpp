#ifndef __FWDPP_EXPERIMENTAL_SAMPLE_DIPLOID_HPP__
#define __FWDPP_EXPERIMENTAL_SAMPLE_DIPLOID_HPP__

namespace KTfwd {
  namespace experimental {
    struct STANDARDMFC
    {
      template<typename T,typename fitness_func>
      inline std::pair<double,std::vector<double> > operator()(const T * diploids,
							       const fitness_func & ff)const
      {
	using diploid_geno_t = typename T::value_type;
	std::vector<double> fitnesses(diploids->size());
	double wbar = 0.;
      
	auto dptr = diploids->begin();
	unsigned N_curr = diploids->size();
	for( unsigned i = 0 ; i < N_curr ; ++i )
	  {
	    (dptr+i)->first->n = 0;
	    (dptr+i)->second->n = 0;
	    fitnesses[i] = fwdpp_internal::diploid_fitness_dispatch(ff,(dptr+i),
								    tags::diploid_type<std::is_base_of<tags::custom_diploid_t,diploid_geno_t>::value>());
	    wbar += fitnesses[i];
	  }
	wbar /= double(diploids->size());
	return std::make_pair(wbar,fitnesses);
      }
    };

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
	      template<typename,typename> class diploid_vector_type,
	      typename mean_fitness_calculator>
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
		   const double & f,
		   const mean_fitness_calculator & mfc)
    {
      assert(N_curr == diploids->size());

      std::for_each( mutations->begin(),mutations->end(),[](typename gamete_type::mutation_type & __m){__m.n=0;});
      //HACK
      auto VV = mfc(diploids,ff);
      double wbar=VV.first;
      auto fitnesses = std::move(VV.second);
      // std::vector<double> fitnesses(diploids->size());
      // double wbar = 0.;
    
      // auto dptr = diploids->begin();
      // for( unsigned i = 0 ; i < N_curr ; ++i )
      //   {
      // 	(dptr+i)->first->n = 0;
      // 	(dptr+i)->second->n = 0;
      // 	fitnesses[i] = fwdpp_internal::diploid_fitness_dispatch(ff,(dptr+i),
      // 								tags::diploid_type<std::is_base_of<tags::custom_diploid_t,diploid_geno_t>::value>());
      // 	wbar += fitnesses[i];
      //   }
      // wbar /= double(diploids->size());

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
	}
      auto dptr = diploids->begin();
      unsigned NREC=0;
      assert(diploids->size()==N_next);
      decltype( gametes->begin() ) p1g1,p1g2,p2g1,p2g2;
      for( unsigned i = 0 ; i < N_next ; ++i )
	{
	  assert(dptr==diploids->begin());
	  assert( (dptr+i) < diploids->end() );
	  size_t p1 = gsl_ran_discrete(r,lookup.get());
	  size_t p2 = (gsl_rng_uniform(r) <= f) ? p1 : gsl_ran_discrete(r,lookup.get());
	  assert(p1<parents.size());
	  assert(p2<parents.size());
	
	  p1g1 = (pptr+typename decltype(pptr)::difference_type(p1))->first;
	  p1g2 = (pptr+typename decltype(pptr)::difference_type(p1))->second;
	  p2g1 = (pptr+typename decltype(pptr)::difference_type(p2))->first;
	  p2g2 = (pptr+typename decltype(pptr)::difference_type(p2))->second;
	
	  NREC += rec_pol(p1g1,p1g2);
	  NREC += rec_pol(p2g1,p2g2);
	
	  (dptr+i)->first = (gsl_rng_uniform(r) <= 0.5) ? p1g1 : p1g2;
	  (dptr+i)->second = (gsl_rng_uniform(r) <= 0.5) ? p2g1 : p2g2;
	
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
#endif
      decltype(gametes->begin()) temp;
      for( auto itr = gametes->begin() ; itr != gametes->end() ;  )
	{
	  if( itr->n == 0 ) //this gamete is extinct and need erasing from the list
	    {
	      temp = itr;
	      ++itr;
	      gametes->erase(temp);
	    }
	  else //gamete remains extant and we adjust mut counts
	    {
	      adjust_mutation_counts(itr,itr->n);
	      ++itr;
	    }
	}
      std::for_each( gametes->begin(),
		     gametes->end(),
		     [&mp](decltype( *(gametes->begin()) ) & __g ) {
		       __g.mutations.erase( std::remove_if(__g.mutations.begin(),__g.mutations.end(),std::cref(mp)),__g.mutations.end() );
		       __g.smutations.erase( std::remove_if(__g.smutations.begin(),__g.smutations.end(),std::cref(mp)),__g.smutations.end() );
		     });
      assert(check_sum(gametes,2*N_next));
      return wbar;
    }
  }
}
