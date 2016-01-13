//  -*- C++ -*-
#ifndef __FWDPP_SAMPLE_DIPLOID_TCC__
#define __FWDPP_SAMPLE_DIPLOID_TCC__

#include <fwdpp/internal/recycling.hpp>
#include <fwdpp/internal/gsl_discrete.hpp>
#include <fwdpp/internal/diploid_fitness_dispatch.hpp>
#include <fwdpp/internal/gamete_lookup_table.hpp>
#include <fwdpp/internal/gamete_cleaner.hpp>
#include <fwdpp/internal/multilocus_rec.hpp>
#include <fwdpp/internal/sample_diploid_helpers.hpp>

namespace KTfwd
{
  //single deme, constant N
  template< typename gamete_type,
	    typename gamete_cont_type_allocator,
	    typename mutation_type,
	    typename mutation_cont_type_allocator,
	    typename diploid_geno_t,
	    typename diploid_vector_type_allocator,
	    typename diploid_fitness_function,
	    typename mutation_model,
	    typename recombination_policy,
	    template<typename,typename> class gamete_cont_type,
	    template<typename,typename> class mutation_cont_type,
	    template<typename,typename> class diploid_vector_type,
	    typename mutation_removal_policy,
	    typename gamete_insertion_policy
	    >
  double
  sample_diploid(gsl_rng * r,
		 gamete_cont_type<gamete_type,gamete_cont_type_allocator > & gametes,
		 diploid_vector_type<diploid_geno_t,diploid_vector_type_allocator> & diploids,
		 mutation_cont_type<mutation_type,mutation_cont_type_allocator > & mutations,
		 std::vector<uint_t> & mcounts,
		 const uint_t & N_curr,
		 const double & mu,
		 const mutation_model & mmodel,
		 const recombination_policy & rec_pol,
		 const diploid_fitness_function & ff,
		 typename gamete_type::mutation_container & neutral,
		 typename gamete_type::mutation_container & selected,
		 const double f,
		 const mutation_removal_policy & mp,
		 const gamete_insertion_policy & gpolicy_mut)
  {
    //run changing N version with N_next == N_curr
    return sample_diploid(r,gametes,diploids,mutations,mcounts,N_curr,N_curr,mu,mmodel,rec_pol,
			  ff,neutral,selected,f,
			  mp,
			  gpolicy_mut);
  }

  //single deme, N changing
  template< typename gamete_type,
	    typename gamete_cont_type_allocator,
	    typename mutation_type,
	    typename mutation_cont_type_allocator,
	    typename diploid_geno_t,
	    typename diploid_vector_type_allocator,
	    typename diploid_fitness_function,
	    typename mutation_model,
	    typename recombination_policy,
	    template<typename,typename> class gamete_cont_type,
	    template<typename,typename> class mutation_cont_type,
	    template<typename,typename> class diploid_vector_type,
	    typename mutation_removal_policy,
	    typename gamete_insertion_policy
	    >
  double
  sample_diploid(gsl_rng * r,
		 gamete_cont_type<gamete_type,gamete_cont_type_allocator > & gametes,
		 diploid_vector_type<diploid_geno_t,diploid_vector_type_allocator> & diploids,
		 mutation_cont_type<mutation_type,mutation_cont_type_allocator > & mutations,
		 std::vector<uint_t> & mcounts,
		 const uint_t & N_curr,
		 const uint_t & N_next,
		 const double & mu,
		 const mutation_model & mmodel,
		 const recombination_policy & rec_pol,
		 const diploid_fitness_function & ff,
		 typename gamete_type::mutation_container & neutral,
		 typename gamete_type::mutation_container & selected,
		 const double f,
		 const mutation_removal_policy & mp,
		 const gamete_insertion_policy & gpolicy_mut)
  {
    //test preconditions in debugging mode
    assert(popdata_sane(diploids,gametes,mcounts));
    assert(mcounts.size()==mutations.size());
    assert(N_curr == diploids.size());
    assert(mcounts.size()==mutations.size());

    std::vector<double> fitnesses(diploids.size());
    double wbar = 0.;
    auto mut_recycling_bin = fwdpp_internal::make_mut_queue(mcounts);
    auto gam_recycling_bin = fwdpp_internal::make_gamete_queue(gametes);
    auto gamete_lookup = fwdpp_internal::gamete_lookup_table(gametes,mutations);
    
    for( uint_t i = 0 ; i < N_curr ; ++i )
      {
	gametes[diploids[i].first].n=gametes[diploids[i].second].n=0;
	fitnesses[i] = fwdpp_internal::diploid_fitness_dispatch(ff,diploids[i],gametes,mutations,
								typename traits::is_custom_diploid_t<diploid_geno_t>::type());
	wbar += fitnesses[i];
      }
    wbar /= double(diploids.size());
#ifndef NDEBUG
    for(const auto & g : gametes) assert(!g.n);
#endif
    fwdpp_internal::gsl_ran_discrete_t_ptr lookup(gsl_ran_discrete_preproc(N_curr,fitnesses.data()));
    const auto parents(diploids); //copy the parents

    //Change the population size
    if( diploids.size() != N_next)
      {
	diploids.resize(N_next);
      }
    assert(diploids.size()==N_next);

    for(auto & dip : diploids)
      {
	size_t p1 = gsl_ran_discrete(r,lookup.get());
	size_t p2 = (f==1. || (f>0. && gsl_rng_uniform(r) < f)) ? p1 : gsl_ran_discrete(r,lookup.get());
	assert(p1<parents.size());
	assert(p2<parents.size());

	size_t p1g1 = parents[p1].first;
	size_t p1g2 = parents[p1].second;
	size_t p2g1 = parents[p2].first;
	size_t p2g2 = parents[p2].second;

	if(gsl_rng_uniform(r)<0.5) std::swap(p1g1,p1g2);
	if(gsl_rng_uniform(r)<0.5) std::swap(p2g1,p2g2);

	dip.first = recombination(gametes,gamete_lookup,gam_recycling_bin,
				  neutral,selected,rec_pol,p1g1,p1g2,mutations).first;
	dip.second = recombination(gametes,gamete_lookup,gam_recycling_bin,
				   neutral,selected,rec_pol,p2g1,p2g2,mutations).first;

	gametes[dip.first].n++;
	gametes[dip.second].n++;

	//now, add new mutations
	dip.first = mutate_gamete_recycle(mut_recycling_bin,gam_recycling_bin,r,mu,gametes,mutations,dip.first,mmodel,gpolicy_mut);
	dip.second = mutate_gamete_recycle(mut_recycling_bin,gam_recycling_bin,r,mu,gametes,mutations,dip.second,mmodel,gpolicy_mut);

	assert( gametes[dip.first].n );
	assert( gametes[dip.second].n );
      }
    assert(check_sum(gametes,2*N_next));
#ifndef NDEBUG
    for(const auto & dip : diploids)
      {
	assert(gametes[dip.first].n>0);
	assert(gametes[dip.first].n<=2*N_next);
	assert(gametes[dip.second].n>0);
	assert(gametes[dip.second].n<=2*N_next);
      }
#endif
    fwdpp_internal::process_gametes(gametes,mutations,mcounts);
    assert(mcounts.size()==mutations.size());
#ifndef NDEBUG
    for(const auto & mc : mcounts)
      {
	assert(mc <= 2*N_next);
      }
#endif
    assert(popdata_sane(diploids,gametes,mcounts));
    fwdpp_internal::gamete_cleaner(gametes,mcounts,2*N_next,typename std::is_same<decltype(mp),KTfwd::remove_nothing >::type());
    return wbar;
  }

  //Metapopulation version of sample_diploid for individual-based simulations and constant N
  template< typename gamete_type,
	    typename mutation_type,
	    typename metapop_diploid_vector_type_allocator,
	    typename gamete_cont_type_allocator,
	    typename mutation_cont_type_allocator,
	    typename diploid_geno_t,
	    typename diploid_vector_type_allocator,
	    typename diploid_fitness_function_container,
	    typename mutation_model,
	    typename recombination_policy,
	    typename migration_policy,
	    template<typename,typename> class gamete_cont_type,
	    template<typename,typename> class mutation_cont_type,
	    template<typename,typename> class diploid_vector_type,
	    template<typename,typename> class metapop_diploid_vector_type,
	    typename mutation_removal_policy,
	    typename gamete_insertion_policy>
  std::vector< double >
  sample_diploid(gsl_rng * r,
		 gamete_cont_type<gamete_type,gamete_cont_type_allocator> & metapop,
		 metapop_diploid_vector_type < diploid_vector_type<diploid_geno_t,diploid_vector_type_allocator>,metapop_diploid_vector_type_allocator > & diploids,
		 mutation_cont_type<mutation_type,mutation_cont_type_allocator > & mutations,
		 std::vector<uint_t> & mcounts,
		 const uint_t * N_curr,
		 const double & mu,
		 const mutation_model & mmodel,
		 const recombination_policy & rec_pol,
		 const diploid_fitness_function_container & ffs,
		 const migration_policy & mig,
		 typename gamete_type::mutation_container & neutral,
		 typename gamete_type::mutation_container & selected,
		 const double * f,
		 const mutation_removal_policy & mp,
		 const gamete_insertion_policy & gpolicy_mut)
  {
    //run changing-N version with no change in N
    return sample_diploid(r,metapop,diploids,mutations,mcounts,N_curr,N_curr,mu,mmodel,rec_pol,
			  ffs,mig,neutral,selected,f,mp,gpolicy_mut);
  }

  //Metapopulation version of sample_diploid for individual-based simulations with changing population size
  template< typename gamete_type,
	    typename mutation_type,
	    typename metapop_diploid_vector_type_allocator,
	    typename gamete_cont_type_allocator,
	    typename mutation_cont_type_allocator,
	    typename diploid_geno_t,
	    typename diploid_vector_type_allocator,
	    typename diploid_fitness_function_container,
	    typename mutation_model,
	    typename recombination_policy,
	    typename migration_policy,
	    template<typename,typename> class gamete_cont_type,
	    template<typename,typename> class mutation_cont_type,
	    template<typename,typename> class diploid_vector_type,
	    template<typename,typename> class metapop_diploid_vector_type,
	    typename mutation_removal_policy,
	    typename gamete_insertion_policy>
  std::vector< double >
  sample_diploid(gsl_rng * r,
		 gamete_cont_type<gamete_type,gamete_cont_type_allocator> & gametes,
		 metapop_diploid_vector_type < diploid_vector_type<diploid_geno_t,diploid_vector_type_allocator>,metapop_diploid_vector_type_allocator > & diploids,
		 mutation_cont_type<mutation_type,mutation_cont_type_allocator > & mutations,
		 std::vector<uint_t> & mcounts,
		 const uint_t * N_curr,
		 const uint_t * N_next,
		 const double & mu,
		 const mutation_model & mmodel,
		 const recombination_policy & rec_pol,
		 const diploid_fitness_function_container & ffs,
		 const migration_policy & mig,
		 typename gamete_type::mutation_container & neutral,
		 typename gamete_type::mutation_container & selected,
		 const double * f,
		 const mutation_removal_policy & mp,
		 const gamete_insertion_policy & gpolicy_mut)
  {
    //get the fitnesses for each diploid in each deme and make the lookup table of parental fitnesses
    using lookup_t = fwdpp_internal::gsl_ran_discrete_t_ptr;
    std::vector<lookup_t> lookups;
    std::vector<double> wbars(diploids.size(),0);
    auto mut_recycling_bin = fwdpp_internal::make_mut_queue(mcounts);
    auto gamete_recycling_bin = fwdpp_internal::make_gamete_queue(gametes);
    auto gamete_lookup = fwdpp_internal::gamete_lookup_table(gametes,mutations);

    //get max N
    auto mN = *std::max_element(N_curr,N_curr+diploids.size());

    std::vector<double> fitnesses(mN);
    std::size_t popi=0;
    for(const auto & dipvec : diploids ) //go over each container of diploids...
      {
	unsigned i=0;
	for(const auto & dip : dipvec) //...and each diploid
	  {
	    fitnesses[i]=fwdpp_internal::diploid_fitness_dispatch(ffs[popi],dip,gametes,mutations,typename traits::is_custom_diploid_t<diploid_geno_t>::type());
	    wbars[popi]+=fitnesses[i];
	    gametes[dip.first].n=gametes[dip.second].n=0;
	    ++i;
	  }
	wbars[popi] /= double(dipvec.size());
	lookups.emplace_back(lookup_t(gsl_ran_discrete_preproc(dipvec.size(),fitnesses.data())));
	++popi;
      }

    assert(lookups.size() == diploids.size());
    //copy diploids into temporary parents
    const auto parents(diploids);

    //Update the diploids, one deme at a time
    for(popi = 0 ; popi < diploids.size() ; ++popi)
      {
	uint_t demesize = *(N_next+popi);
	if(demesize != *(N_curr+popi))
	  {
	    diploids[popi].resize(demesize);
	  }
	for(auto & dip : diploids[popi])
	  {
	    /* Figure out if parent 1 is migrant or not.
	       
	       A migration policy takes the current deme (popindex) as
	       an argument.  It returns popindex if there is no migration,
	       else it returns the index of the deme of a migrant parent
	    */
	    std::size_t deme_p1 = mig(popi),deme_p2=popi;

	    //Figure out who the parents are
	    std::size_t p1 = gsl_ran_discrete(r,lookups[deme_p1].get()),p2;

	    /*
	      If the individual is not inbred, then we pick a
	      deme from the migration policy for parent 2
	    */
	    if( f != nullptr && ( *(f + popi)==1. || (*(f + popi)>0. && gsl_rng_uniform(r) < *(f + popi)) ) ) //individual is inbred
	      {
		p2=p1;
	      }
	    else
	      {
		//apply migration policy to figure out parental deme for parent #2
		deme_p2 = mig(popi);
		p2 = gsl_ran_discrete(r,lookups[deme_p2].get());
	      }

	    std::size_t p1g1 = parents[deme_p1][p1].first;
	    std::size_t p1g2 = parents[deme_p1][p1].second;
	    std::size_t p2g1 = parents[deme_p2][p2].first;
	    std::size_t p2g2 = parents[deme_p2][p2].second;

	    if(gsl_rng_uniform(r)<0.5)std::swap(p1g1,p1g2);
	    if(gsl_rng_uniform(r)<0.5)std::swap(p2g1,p2g2);

	    dip.first = recombination(gametes,gamete_lookup,gamete_recycling_bin,
				      neutral,selected,rec_pol,p1g1,p1g2,mutations).first;
	    dip.second = recombination(gametes,gamete_lookup,gamete_recycling_bin,
				       neutral,selected,rec_pol,p2g1,p2g2,mutations).first;
	    
	    gametes[dip.first].n++;
	    gametes[dip.second].n++;
	    
	    //now, add new mutations
	    dip.first = mutate_gamete_recycle(mut_recycling_bin,gamete_recycling_bin,r,mu,gametes,mutations,dip.first,mmodel,gpolicy_mut);
	    dip.second = mutate_gamete_recycle(mut_recycling_bin,gamete_recycling_bin,r,mu,gametes,mutations,dip.second,mmodel,gpolicy_mut);
	  }
      }
    fwdpp_internal::process_gametes(gametes,mutations,mcounts);
    assert(mcounts.size()==mutations.size());
    fwdpp_internal::gamete_cleaner(gametes,mcounts,
				   2*std::accumulate(N_next,N_next+diploids.size(),uint_t(0)),
				   typename std::is_same<decltype(mp),KTfwd::remove_nothing >::type());
    assert( check_sum(gametes,2*std::accumulate(N_next,N_next+diploids.size(),uint_t(0))) );
    return wbars;
  }
  
  //Multi-locus API
  //single deme, N changing
  template< typename diploid_geno_t,
	    typename gamete_type,
	    typename gamete_cont_type_allocator,
	    typename mutation_type,
	    typename mutation_cont_type_allocator,
	    typename diploid_vector_type_allocator,
	    typename locus_vector_type_allocator,
	    typename diploid_fitness_function,
	    typename mutation_model_container,
	    typename recombination_policy_container,
	    typename bw_locus_rec_fxn,
	    template<typename,typename> class gamete_cont_type,
	    template<typename,typename> class mutation_cont_type,
	    template<typename,typename> class diploid_vector_type,
	    template<typename,typename> class locus_vector_type,
	    typename mutation_removal_policy,
	    typename gamete_insertion_policy
	    >
  double
  sample_diploid(gsl_rng * r,
		 gamete_cont_type<gamete_type, gamete_cont_type_allocator> & gametes,
		 diploid_vector_type<locus_vector_type< diploid_geno_t, locus_vector_type_allocator>, diploid_vector_type_allocator> & diploids,
		 mutation_cont_type<mutation_type,mutation_cont_type_allocator > & mutations,
		 std::vector<uint_t> & mcounts,
		 const uint_t & N_curr,
		 const uint_t & N_next,
		 const double * mu,
		 const mutation_model_container & mmodel,
		 const recombination_policy_container & rec_policies,
		 const double * r_between_loci,
		 const bw_locus_rec_fxn & blrf,
		 const diploid_fitness_function & ff,
		 typename gamete_type::mutation_container & neutral,
		 typename gamete_type::mutation_container & selected,
		 const double & f,
		 const mutation_removal_policy & mp,
		 const gamete_insertion_policy & gpolicy_mut)
  {
    assert(mcounts.size()==mutations.size());
#ifndef NDEBUG
    assert(diploids.size()==N_curr);
#endif
    //Vector of parental fitnesses
    std::vector<double> fitnesses(N_curr);
    double wbar = 0.;
    auto mut_recycling_bin = fwdpp_internal::make_mut_queue(mcounts);
    auto gamete_recycling_bin = fwdpp_internal::make_gamete_queue(gametes);
    auto gamete_lookup = fwdpp_internal::gamete_lookup_table(gametes,mutations);
    //Go over parents
    for( uint_t i = 0 ; i < N_curr ; ++i )
      {
	//set parental gamete counts to 0 for each locus
	for( auto locus : diploids[i])
	  {
	    gametes[locus.first].n=gametes[locus.second].n=0;
	  }
	//Calculate the fitness of this parent
	fitnesses[i] += ff( diploids[i] );
	//increment pop. mean fitness
	wbar += fitnesses[i];
      }
    wbar /= double(diploids.size());
#ifndef NDEBUG
    /*
      If we are debugging, let's make sure that every gamete has been set to n = 0.
      Rationale for check:  if we are failing to update data types properly, then
      it is possible that the "gamete pool" contains items not carried by any diploids.
      If so, this assertion will fail.
    */
    for(const auto & g : gametes) assert(!g.n);
#endif

    fwdpp_internal::gsl_ran_discrete_t_ptr lookup(gsl_ran_discrete_preproc(fitnesses.size(),fitnesses.data()));

    const auto parents(diploids); //Copy the parents.  Exact copy of diploids--same fitnesses, etc.

    //Change the population size and reset dptr to avoid iterator invalidation
    if( diploids.size() != N_next)
      {
	diploids.resize(N_next);
      }

    assert(diploids.size()==N_next);

    for(auto & dip : diploids)
      {
	std::size_t p1 = gsl_ran_discrete(r,lookup.get());
	assert(p1<N_curr);
	std::size_t p2 = (f==1.||(f>0. && gsl_rng_uniform(r)<f)) ? p1 : gsl_ran_discrete(r,lookup.get());
	assert(p2<N_curr);
	dip=fwdpp_internal::multilocus_rec_mut(r,parents[p1],parents[p2],
					       mut_recycling_bin,gamete_recycling_bin,gamete_lookup,
					       rec_policies,blrf,r_between_loci,
					       ((gsl_rng_uniform(r)<0.5)?1:0),
					       ((gsl_rng_uniform(r)<0.5)?1:0),
					       gametes,mutations,neutral,selected,mu,mmodel,gpolicy_mut
					       );
	
      }
    fwdpp_internal::process_gametes(gametes,mutations,mcounts);
    fwdpp_internal::gamete_cleaner(gametes,mcounts,2*N_next,
				   typename std::is_same<decltype(mp),KTfwd::remove_nothing >::type());
    return wbar;
  }

  //single deme, constant N
  template< typename diploid_geno_t,
	    typename gamete_type,
	    typename gamete_cont_type_allocator,
	    typename mutation_type,
	    typename mutation_cont_type_allocator,
	    typename diploid_vector_type_allocator,
	    typename locus_vector_type_allocator,
	    typename diploid_fitness_function,
	    typename mutation_model_container,
	    typename recombination_policy_container,
	    typename bw_locus_rec_fxn,
	    template<typename,typename> class gamete_cont_type,
	    template<typename,typename> class mutation_cont_type,
	    template<typename,typename> class diploid_vector_type,
	    template<typename,typename> class locus_vector_type,
	    typename mutation_removal_policy,
	    typename gamete_insertion_policy>
  double
  sample_diploid(gsl_rng * r,
		      gamete_cont_type<gamete_type, gamete_cont_type_allocator> & gametes,
		      diploid_vector_type<locus_vector_type< diploid_geno_t, locus_vector_type_allocator>, diploid_vector_type_allocator> & diploids,
		      mutation_cont_type<mutation_type,mutation_cont_type_allocator > & mutations,
		      std::vector<uint_t> & mcounts,
		      const uint_t & N,
		      const double * mu,
		      const mutation_model_container & mmodel,
		      const recombination_policy_container & rec_policies,
		      const double * r_between_loci,
		      const bw_locus_rec_fxn & blrf,
		      const diploid_fitness_function & ff,
		      typename gamete_type::mutation_container & neutral,
		      typename gamete_type::mutation_container & selected,
		      const double & f,
		      const mutation_removal_policy & mp,
		      const gamete_insertion_policy & gpolicy_mut)
  {
    return sample_diploid(r,gametes,diploids,mutations,mcounts,N,N,mu,mmodel,rec_policies,r_between_loci,blrf,ff,neutral,selected,f,mp,gpolicy_mut);
  }
}

#endif

