#ifndef __FWDPP_EXPERIMENTAL_SAMPLE_DIPLOID_HPP__
#define __FWDPP_EXPERIMENTAL_SAMPLE_DIPLOID_HPP__

namespace KTfwd {
  namespace experimental {
    /*!
      \brief Abstraction of the standard Wright-Fisher sampling process
     */
    struct standardWFrules
    {
      mutable double wbar;
      mutable std::vector<double> fitnesses;

      mutable fwdpp_internal::gsl_ran_discrete_t_ptr lookup;
      //! \brief Constructor
      standardWFrules() : wbar(0.),fitnesses(std::vector<double>()),lookup(fwdpp_internal::gsl_ran_discrete_t_ptr(nullptr))
      {
      }

      //! \brief The "fitness manager"
      template<typename T,typename fitness_func>
      void w(const T * diploids,
	     const fitness_func & ff)const
      {
	using diploid_geno_t = typename T::value_type;
	unsigned N_curr = diploids->size();
	if(fitnesses.size() < N_curr) fitnesses.resize(N_curr);
	wbar = 0.;
      
	auto dptr = diploids->begin();

	for( unsigned i = 0 ; i < N_curr ; ++i )
	  {
	    (dptr+i)->first->n = 0;
	    (dptr+i)->second->n = 0;
	    fitnesses[i] = fwdpp_internal::diploid_fitness_dispatch(ff,(dptr+i),
								    tags::diploid_type<std::is_base_of<tags::custom_diploid_t,diploid_geno_t>::value>());
	    wbar += fitnesses[i];
	  }
	wbar /= double(diploids->size());

	/*!
	  Black magic alert:
	  fwdpp_internal::gsl_ran_discrete_t_ptr contains a std::unique_ptr wrapping the GSL pointer.
	  This type has its own deleter, which is convenient, because
	  operator= for unique_ptrs automagically calls the deleter before assignment!
	  Details: http://www.cplusplus.com/reference/memory/unique_ptr/operator=
	*/
	lookup = fwdpp_internal::gsl_ran_discrete_t_ptr(gsl_ran_discrete_preproc(N_curr,&fitnesses[0]));
      }

      //! \brief Pick parent one
      inline size_t pick1(gsl_rng * r) const
      {
	return gsl_ran_discrete(r,lookup.get());
      }

      //! \brief Pick parent 2.  Parent 1's data are passed along for models where that is relevant
      template<typename diploid_itr_t>
      inline size_t pick2(gsl_rng * r, const size_t & p1, diploid_itr_t p1_itr, const double & f ) const
      {
	//static asserts suppress hideously-long compiler warnings on GCC
	static_assert( std::is_const< typename std::remove_pointer<typename decltype(p1_itr)::pointer>::type >::value , "p1_itr must point to const data");
	return (gsl_rng_uniform(r) <= f) ? p1 : gsl_ran_discrete(r,lookup.get());
      }

      //! \brief Update some property of the offspring based on properties of the parents
      template<typename offspring_itr_t, typename parent_itr_t>
      void update(gsl_rng * r, offspring_itr_t offspring, parent_itr_t p1_itr, parent_itr_t p2_itr ) const
      {
	//static asserts suppress hideously-long compiler warnings on GCC
	static_assert( std::is_pointer<decltype(r)>::value, "r must be a pointer" );
	static_assert(!std::is_const< typename std::remove_pointer<typename decltype(offspring)::pointer>::type >::value , "offpsring must point to non-const data");
	static_assert( std::is_const< typename std::remove_pointer<typename decltype(p1_itr)::pointer>::type >::value , "p1_itr must point to const data");
	static_assert( std::is_const< typename std::remove_pointer<typename decltype(p2_itr)::pointer>::type >::value , "p2_itr must point to const data");
	return;
      }

    };

    //single deme, N changing

    //! \brief Experimental variant where the population rules are implemented via an external policy
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
	      typename popmodel_rules = standardWFrules>
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
		   const double & f = 0.,
		   const popmodel_rules & pmr = popmodel_rules())
    {
      assert(N_curr == diploids->size());

      std::for_each( mutations->begin(),mutations->end(),[](typename gamete_type::mutation_type & __m){__m.n=0;});

      pmr.w(diploids,ff);

#ifndef NDEBUG
      std::for_each(gametes->cbegin(),gametes->cend(),[](decltype((*gametes->cbegin())) __g) {
	  assert( !__g.n ); } );
#endif
      auto parents(*diploids); //copy the parents
      auto pptr = parents.cbegin();
    
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
	  //Pick parent 1
	  size_t p1 = pmr.pick1(r);
	  //Pick parent 2
	  size_t p2 = pmr.pick2(r,p1,pptr+typename decltype(pptr)::difference_type(p1),f);
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

	  pmr.update(r,(dptr+i),pptr+typename decltype(pptr)::difference_type(p1),pptr+typename decltype(pptr)::difference_type(p2));
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
      return pmr.wbar;
    }

    //single deme, N constant

    //! \brief Experimental variant where the population rules are implemented via an external policy
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
	      typename popmodel_rules = standardWFrules>
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
		   const double & f = 0.,
		   const popmodel_rules & pmr = popmodel_rules())
    {
      return experimental::sample_diploid(r,gametes,diploids,mutations,N_curr,N_curr,mu,mmodel,rec_pol,mpolicy,gpolicy_mut,ff,mp,f,pmr);
    }
  }
}
#endif
