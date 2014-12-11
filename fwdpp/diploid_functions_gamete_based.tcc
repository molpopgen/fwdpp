//  -*- C++ -*- 
#ifndef __DIPLOID_FUNCTIONS_GAMETE_BASED_HPP__
#define __DIPLOID_FUNCTIONS_GAMETE_BASED_HPP__

#if defined(HAVE_BOOST_VECTOR) && !defined(USE_STANDARD_CONTAINERS)
#include <boost/container/vector.hpp>
#endif

#include <fwdpp/internal/recombination_common.hpp>
#include <fwdpp/internal/gsl_discrete.hpp>
#include <vector>

namespace KTfwd
{
  template< typename gamete_type,
	    typename diploid_fitness_function,
	    typename floating_type,
	    typename vector_type_allocator,
	    template<typename,typename> class vector_type>
  inline void calc_mean_fitness( floating_type *wbar,
				 floating_type * expected_gamete_frequencies,
				 const vector_type<gamete_type,vector_type_allocator > * gametes, 
				 const unsigned & twoN, const diploid_fitness_function & ff,
				 const double & f)
  {
    assert(!gametes->empty());
    typename vector_type<gamete_type,vector_type_allocator >::const_iterator ibeg=gametes->begin(),jbeg;
    double pqw,p,q,w;
    unsigned j;
    for(unsigned i=0;i<gametes->size();++i)
      {
	if(ibeg->n)
	  {
	    p = double(ibeg->n)/double(twoN);
	    if ( (w= ff( ibeg, ibeg )) >0. )
	      {
		pqw = (p*p*(1.-f) + p*f)*w;
		*wbar += pqw;
		*(expected_gamete_frequencies+i) += pqw;
	      }
	  
	    j=i+1;
	    jbeg=ibeg+1;
	    for( ; j < gametes->size() ; ++j)
	      {
		if(jbeg->n)
		  {
		    q = double(jbeg->n)/double(twoN);
		    if((w = ff( ibeg, jbeg ))>0.)
		      {
			pqw=p*q*(1.-f)*w;
			*wbar += 2.*pqw;
			*(expected_gamete_frequencies+i) += pqw;
			*(expected_gamete_frequencies+j) += pqw;
		      }
		  }
		++jbeg;
	      }
	  }
	++ibeg;
      }
  }

  //SINGLE-POP VERSIONS
  template< typename gamete_type,
	    typename vector_type_allocator,
	    template<typename,typename> class vector_type,
	    typename diploid_fitness_function,
	    typename mutation_removal_policy>
  double sample_diploid(gsl_rng * r, vector_type<gamete_type,vector_type_allocator > * gametes, 
			const unsigned & twoN, 
			const diploid_fitness_function & ff,
			const mutation_removal_policy & mp,
			const double & f)
	    {
              static_assert( std::is_base_of<mutation_base,typename gamete_type::mutation_type>::value,
                             "gamete_type::mutation_value must inherit from KTfwd::mutation_base" );
	      typedef gamete_base< typename gamete_type::mutation_type, typename gamete_type::mutation_list_type > gamete_base_type;
	      static_assert( std::is_base_of<gamete_base_type,gamete_type>::value ||
                             std::is_same<gamete_base_type,gamete_type>::value,
                             "gamete_type must be, or inherit from, KTfwd::gamete_base<mutation_type,mutation_list_type>" );
	      double wbar = 0.;
	      std::vector<double> expected_gamete_frequencies(gametes->size(),0.);
	      calc_mean_fitness(&wbar,&expected_gamete_frequencies[0],gametes,twoN,ff,f);
	      //normalize expected frequqncies
	      std::transform(expected_gamete_frequencies.begin(),expected_gamete_frequencies.end(),
			     expected_gamete_frequencies.begin(),std::bind2nd(std::divides<double>(),wbar));
	      
	      multinomial_sample_gametes(r, gametes->begin(),gametes->size(),&expected_gamete_frequencies[0],twoN);
    
	      update_gamete_list(gametes,twoN,mp);
	      return wbar;
	    }

template< typename gamete_type,
	  typename vector_type_allocator,
	  template<typename,typename> class vector_type,
	  typename diploid_fitness_function,
	  typename mutation_removal_policy>
double sample_diploid(gsl_rng * r, vector_type<gamete_type,vector_type_allocator > * gametes, 
		      const unsigned & twoN_curr, const unsigned & twoN_next,
		      const diploid_fitness_function & ff,
		      const mutation_removal_policy & mp,
		      const double & f)
//hermaphroditic diploids
{
  static_assert( std::is_base_of<mutation_base,typename gamete_type::mutation_type>::value,
                 "gamete_type::mutation_value must inherit from KTfwd::mutation_base" );
  typedef gamete_base< typename gamete_type::mutation_type, typename gamete_type::mutation_list_type > gamete_base_type;
  static_assert( std::is_base_of<gamete_base_type,gamete_type>::value ||
                 std::is_same<gamete_base_type,gamete_type>::value,
                 "gamete_type must be, or inherit from, KTfwd::gamete_base<mutation_type,mutation_list_type>" );
  double wbar = 0.;
  std::vector<double> expected_gamete_frequencies(gametes->size(),0.);

  calc_mean_fitness(&wbar,&expected_gamete_frequencies[0],gametes,twoN_curr,ff,f);
  //normalize expected frequqncies
  std::transform(expected_gamete_frequencies.begin(),expected_gamete_frequencies.end(),
		 expected_gamete_frequencies.begin(),std::bind2nd(std::divides<double>(),wbar));

  multinomial_sample_gametes(r, gametes->begin(),gametes->size(),&expected_gamete_frequencies[0],twoN_next);
  update_gamete_list(gametes,twoN_next,mp);
  return wbar;
}

  //METAPOP VERSIONS
  template< typename gamete_type,
	    typename vector_type_allocator1,
	    typename vector_type_allocator2,
	    template<typename,typename> class vector_type1,
	    template<typename,typename> class vector_type2,
	    typename diploid_fitness_function_container,
	    typename mutation_removal_policy>
  std::vector<double> sample_diploid(gsl_rng * r, 
				     vector_type2< vector_type1<gamete_type, vector_type_allocator1 >, vector_type_allocator2 > * metapop,
				     const unsigned * twoNs,
				     const unsigned & metapopsize,
				     const diploid_fitness_function_container & ffs,
				     const mutation_removal_policy & mp)
  //hermaphroditic diploids
	    {
	      static_assert( std::is_base_of<mutation_base,typename gamete_type::mutation_type>::value,
                             "gamete_type::mutation_value must inherit from KTfwd::mutation_base" );
	      typedef gamete_base< typename gamete_type::mutation_type, typename gamete_type::mutation_list_type > gamete_base_type;
              static_assert( std::is_base_of<gamete_base_type,gamete_type>::value ||
                             std::is_same<gamete_base_type,gamete_type>::value,
                             "gamete_type must be, or inherit from, KTfwd::gamete_base<mutation_type,mutation_list_type>" );
	      assert( metapop->size() == ffs.size() );
	      assert(std::accumulate(twoNs,twoNs + metapop->size(),0u) == metapopsize);
	      std::vector<double> wbars;

	      typedef typename vector_type2< vector_type1<gamete_type, vector_type_allocator1 >, vector_type_allocator2 >::iterator mpop_itr;

	      unsigned i = 0;
	      for(mpop_itr pop = metapop->begin() ; pop < metapop->end() ; ++pop,++i)
		{
		  double wbar = 0.;
		  std::vector<double> expected_gamete_frequencies(pop->size(),0.);
		  calc_mean_fitness(&wbar,&expected_gamete_frequencies[0],
				    (&*pop),
				    *(twoNs+i),ffs[i],0.);
		  wbars.push_back(wbar);
		  std::transform(expected_gamete_frequencies.begin(),expected_gamete_frequencies.end(),
				 expected_gamete_frequencies.begin(),std::bind2nd(std::divides<double>(),wbar));
		  multinomial_sample_gametes(r, pop->begin(),pop->size(),&expected_gamete_frequencies[0],*(twoNs+i));
		}
	      
	      for(mpop_itr pop = metapop->begin() ; pop < metapop->end() ; ++pop)
		{
		  update_gamete_list(&*pop,metapopsize,mp);
		}
	      return wbars;
	    }

template< typename gamete_type,
	  typename vector_type_allocator1,
	  typename vector_type_allocator2,
	  template<typename,typename> class vector_type1,
	  template<typename,typename> class vector_type2,
	  typename diploid_fitness_function_container,
	  typename mutation_removal_policy>
std::vector<double> sample_diploid(gsl_rng * r, 
				   vector_type2< vector_type1<gamete_type, vector_type_allocator1 >, vector_type_allocator2 > * metapop,
				   const unsigned * twoNs,
				   const unsigned * twoNs_next,
				   const unsigned & metapopsize,
				   const diploid_fitness_function_container & ffs,
				   const mutation_removal_policy & mp)
//hermaphroditic diploids
{
  static_assert( std::is_base_of<mutation_base,typename gamete_type::mutation_type>::value,
                 "gamete_type::mutation_value must inherit from KTfwd::mutation_base" );
  typedef gamete_base< typename gamete_type::mutation_type, typename gamete_type::mutation_list_type > gamete_base_type;
  static_assert( std::is_base_of<gamete_base_type,gamete_type>::value ||
		std::is_same<gamete_base_type,gamete_type>::value,
		"gamete_type must be, or inherit from, KTfwd::gamete_base<mutation_type,mutation_list_type>" );
  assert( metapop->size() == ffs.size() );
  assert(std::accumulate(twoNs,twoNs + metapop->size(),0) == metapopsize);
  std::vector<double> wbars;

  typedef typename vector_type2< vector_type1<gamete_type, vector_type_allocator1 >, vector_type_allocator2 >::iterator mpop_itr;

  unsigned i = 0;
  for(mpop_itr pop = metapop->begin() ; pop < metapop->end() ; ++pop,++i)
    {
      double wbar = 0.;
      std::vector<double> expected_gamete_frequencies(pop->size(),0.);
      calc_mean_fitness(&wbar,&expected_gamete_frequencies[0],
			(&*pop),
			*(twoNs+i),ffs[i],0.);
      wbars.push_back(wbar);
      std::transform(expected_gamete_frequencies.begin(),expected_gamete_frequencies.end(),
		     expected_gamete_frequencies.begin(),std::bind2nd(std::divides<double>(),wbar));
      multinomial_sample_gametes(r, pop->begin(),pop->size(),&expected_gamete_frequencies[0],*(twoNs_next+i));
    }
	      
  for(mpop_itr pop = metapop->begin() ; pop < metapop->end() ; ++pop)
    {
      update_gamete_list(&*pop,std::accumulate(twoNs_next,twoNs_next+metapop->size(),0u),mp);
    }
  return wbars;
}

  template< typename gamete_type,
	    typename vector_type_allocator1,
	    typename vector_type_allocator2,
	    template<typename,typename> class vector_type1,
	    template<typename,typename> class vector_type2,
	    typename diploid_fitness_function_container,
	    typename mutation_removal_policy>
  std::vector<double> sample_diploid(gsl_rng * r, 
				     vector_type2< vector_type1<gamete_type, vector_type_allocator1 >, vector_type_allocator2 > * metapop,
				     const unsigned * twoNs,
				     const unsigned & metapopsize,
				     const diploid_fitness_function_container & ffs,
				     const mutation_removal_policy & mp,
				     const double * fs)
  //hermaphroditic diploids
	    {
	      static_assert( std::is_base_of<mutation_base,typename gamete_type::mutation_type>::value,
			     "gamete_type::mutation_value must inherit from KTfwd::mutation_base" );
	      typedef gamete_base< typename gamete_type::mutation_type, typename gamete_type::mutation_list_type > gamete_base_type;
	      static_assert( std::is_base_of<gamete_base_type,gamete_type>::value ||
                             std::is_same<gamete_base_type,gamete_type>::value,
                             "gamete_type must be, or inherit from, KTfwd::gamete_base<mutation_type,mutation_list_type>" );
	      assert( metapop->size() == ffs.size() );
	      assert(std::accumulate(twoNs,twoNs + metapop->size(),0) == metapopsize);
	      std::vector<double> wbars;

	      typedef typename vector_type2< vector_type1<gamete_type, vector_type_allocator1 >, vector_type_allocator2 >::iterator mpop_itr;

	      unsigned i = 0;
	      for(mpop_itr pop = metapop->begin() ; pop < metapop->end() ; ++pop,++i)
		{
		  double wbar = 0.;
		  std::vector<double> expected_gamete_frequencies(pop->size(),0.);
		  calc_mean_fitness(&wbar,&expected_gamete_frequencies[0],
				    (&*pop),
				    *(twoNs+i),ffs[i],*(fs+i));
		  wbars.push_back(wbar);
		  std::transform(expected_gamete_frequencies.begin(),expected_gamete_frequencies.end(),
				 expected_gamete_frequencies.begin(),std::bind2nd(std::divides<double>(),wbar));
		  multinomial_sample_gametes(r, pop->begin(),pop->size(),&expected_gamete_frequencies[0],*(twoNs+i));
		}
	      
	      for(mpop_itr pop = metapop->begin() ; pop < metapop->end() ; ++pop)
		{
		  update_gamete_list(&*pop,metapopsize,mp);
		}
	      return wbars;
	    }

template< typename gamete_type,
	  typename vector_type_allocator1,
	  typename vector_type_allocator2,
	  template<typename,typename> class vector_type1,
	  template<typename,typename> class vector_type2,
	  typename diploid_fitness_function_container,
	  typename mutation_removal_policy>
std::vector<double> sample_diploid(gsl_rng * r, 
				   vector_type2< vector_type1<gamete_type, vector_type_allocator1 >, vector_type_allocator2 > * metapop,
				   const unsigned * twoNs,
				   const unsigned * twoNs_next,
				   const unsigned & metapopsize,
				   const diploid_fitness_function_container & ffs,
				   const mutation_removal_policy & mp,
				   const double * fs)
//hermaphroditic diploids
{
  static_assert( std::is_base_of<mutation_base,typename gamete_type::mutation_type>::value,
                 "gamete_type::mutation_value must inherit from KTfwd::mutation_base" );
  typedef gamete_base< typename gamete_type::mutation_type, typename gamete_type::mutation_list_type > gamete_base_type;
  static_assert( std::is_base_of<gamete_base_type,gamete_type>::value ||
                 std::is_same<gamete_base_type,gamete_type>::value,
                 "gamete_type must be, or inherit from, KTfwd::gamete_base<mutation_type,mutation_list_type>" );
  assert( metapop->size() == ffs.size() );
  assert(std::accumulate(twoNs,twoNs + metapop->size(),0) == metapopsize);
  std::vector<double> wbars;

  typedef typename vector_type2< vector_type1<gamete_type, vector_type_allocator1 >, vector_type_allocator2 >::iterator mpop_itr;

  unsigned i = 0;
  for(mpop_itr pop = metapop->begin() ; pop < metapop->end() ; ++pop,++i)
    {
      double wbar = 0.;
      std::vector<double> expected_gamete_frequencies(pop->size(),0.);
      calc_mean_fitness(&wbar,&expected_gamete_frequencies[0],
			(&*pop),
			*(twoNs+i),ffs[i],*(fs+i));
      wbars.push_back(wbar);
      std::transform(expected_gamete_frequencies.begin(),expected_gamete_frequencies.end(),
		     expected_gamete_frequencies.begin(),std::bind2nd(std::divides<double>(),wbar));
      multinomial_sample_gametes(r, pop->begin(),pop->size(),&expected_gamete_frequencies[0],*(twoNs_next+i));
    }
	      
  for(mpop_itr pop = metapop->begin() ; pop < metapop->end() ; ++pop)
    {
      update_gamete_list(&*pop,std::accumulate(twoNs_next,twoNs_next+metapop->size(),0u));
    }
  return wbars;
}

template< typename gamete_type,
	  typename vector_type_allocator,
	  template<typename,typename> class vector_type,
	  typename map_function>
unsigned recombine(gsl_rng * r, vector_type<gamete_type,vector_type_allocator > * gametes,
		   const unsigned & twoN, const double & littler,const map_function & mf,
		   const double & f)
{
  static_assert( std::is_base_of<mutation_base,typename gamete_type::mutation_type>::value,
                 "gamete_type::mutation_value must inherit from KTfwd::mutation_base" );
  typedef gamete_base< typename gamete_type::mutation_type, typename gamete_type::mutation_list_type > gamete_base_type;
  static_assert( std::is_base_of<gamete_base_type,gamete_type>::value ||
                 std::is_same<gamete_base_type,gamete_type>::value,
                 "gamete_type must be, or inherit from, KTfwd::gamete_base<mutation_type,mutation_list_type>" );
  typedef typename gamete_type::mutation_container mcont;
  typedef typename vector_type<gamete_type, vector_type_allocator>::iterator vtype_iterator;
  //"individual-based recombination"
  //1. determine # recombinants in whole pop
  //2. assign them randomly to pairs of gametes chosen from the gametes list
  //3. if they are the same gamete, ignore, else do recombination thing.

  //calc freq of homozygotes
  
  double fAA=0.;
#if defined(HAVE_BOOST_VECTOR) && !defined(USE_STANDARD_CONTAINERS)
  boost::container::vector<double> gcounts(gametes->size(),0);
#else
  std::vector<double> gcounts(gametes->size(),0);
#endif
  const double twoNsq = double(twoN)*double(twoN);
  for(vtype_iterator i=gametes->begin() ;i<gametes->end();++i)
    {
      gcounts[decltype(gcounts)::size_type(i - gametes->begin())] = i->n;
      fAA += ( (double(i->n)*double(i->n)*(1.-f))/twoNsq + (double(i->n)/twoN)*f );
    }
  assert(gcounts.size()==gametes->size());
  using lookup_t = fwdpp_internal::gsl_ran_discrete_t_ptr;
  lookup_t lookup(gsl_ran_discrete_preproc(gcounts.size(),&gcounts[0]));
  unsigned NRECS = unsigned(gsl_ran_poisson(r,(1.-fAA)*double(twoN)*littler));
#ifndef NDEBUG
  unsigned doublecheck=NRECS;
  unsigned NRECS_DONE = 0;
  decltype(gametes->size()) ncurrent_classes = gametes->size();
#endif
  vtype_iterator ibeg,jbeg;
  unsigned NEXTINCT=0;
  size_t NRECS_i;

  while(NRECS > 0)
    {
      size_t ith = gsl_ran_discrete(r,lookup.get());
      while(gcounts[ith]==0)
	{
	  ith = gsl_ran_discrete(r,lookup.get());
	}
      size_t jth = gsl_ran_discrete(r,lookup.get());
      while(ith==jth || gcounts[jth]==0)
	{
	  jth = gsl_ran_discrete(r,lookup.get());
	}
      ibeg = (gametes->begin()+typename std::iterator_traits<decltype(gametes->begin())>::difference_type(ith));
      jbeg = (gametes->begin()+typename std::iterator_traits<decltype(gametes->begin())>::difference_type(jth));
      assert( ibeg != gametes->end() &&
	      jbeg != gametes->end() );
      assert(ibeg->n >= 1);
      assert(jbeg->n >= 1);
      decltype(ibeg->mutations.size()) nm1=ibeg->mutations.size()+ibeg->smutations.size(),
	nm2=jbeg->mutations.size()+jbeg->smutations.size();
      //if one gametes carries 0 mutations, and the other carries 1,
      //it is impossible for recombination to generate any new types,
      //so we skip those cases...
      if(!(std::min(nm1,nm2)==0 && std::max(nm1,nm2)==1))
	{
	  //initialize a discrete distribution from 1 to nrecs that is a truncated Poisson
#if defined(HAVE_BOOST_VECTOR) && !defined(USE_STANDARD_CONTAINERS)
	  boost::container::vector<double> breakpoint_dist;
#else
	  std::vector<double> breakpoint_dist;
#endif
	  for(unsigned i = 1 ; i <= NRECS ; ++i)
	    {
	      breakpoint_dist.push_back( gsl_ran_poisson_pdf(i,2.*littler) );
	    }
	  lookup_t BPLOOKUP( gsl_ran_discrete_preproc (size_t(NRECS), &breakpoint_dist[0]));
	  NRECS_i = gsl_ran_discrete(r,BPLOOKUP.get())+1;
	  //assign breakpoints
#if defined(HAVE_BOOST_VECTOR) && !defined(USE_STANDARD_CONTAINERS)
	  boost::container::vector<double> pos;
#else
	  std::vector<double> pos;
#endif
	  for(unsigned i = 0 ; i < NRECS_i ; ++i)
	    {
	      pos.push_back(mf());
	    }
	  std::sort(pos.begin(),pos.end());
	  pos.push_back(std::numeric_limits<double>::max());
	  ibeg->n--;
	  jbeg->n--;
	  gcounts[ith]--;
	  gcounts[jth]--;
	  NEXTINCT += (ibeg->n==0) ? 1u : 0u;
	  NEXTINCT += (jbeg->n==0) ? 1u : 0u;
	  gamete_type new_gamete1(1, mcont(),mcont()),new_gamete2(new_gamete1);
	  new_gamete1.mutations.reserve(ibeg->mutations.size()+jbeg->mutations.size());
	  new_gamete1.smutations.reserve(ibeg->smutations.size()+jbeg->smutations.size());
	  new_gamete2.mutations.reserve(ibeg->mutations.size()+jbeg->mutations.size());
	  new_gamete2.smutations.reserve(ibeg->smutations.size()+jbeg->smutations.size());

	  fwdpp_internal::recombine_gametes(pos,ibeg,jbeg,new_gamete1,new_gamete2);

	  vtype_iterator newpos = update_if_exists_insert(new_gamete1,gametes);
	  assert(newpos->n <= twoN);
	  auto dist = decltype(gcounts.size())(newpos - gametes->begin());
	  if( dist < gcounts.size() )
	    {
	      gcounts[dist]++;
	      assert( gcounts[dist] <= twoN );
	    }
	  newpos = update_if_exists_insert(new_gamete2,gametes);
	  assert(newpos->n <= twoN);
	  dist = decltype(dist)(newpos - gametes->begin());
	  if( dist < gcounts.size() )
	    {
	      gcounts[dist]++;
	      assert( gcounts[dist] <= twoN );
	    }
#ifndef NDEBUG
	  unsigned dummy=0;
	  for( vtype_iterator test = gametes->begin();
	       test!=gametes->begin()+typename std::iterator_traits<decltype(gametes->begin())>::difference_type(ncurrent_classes) ;
	       ++test,++dummy )
	    {
	      assert(test->n == gcounts[dummy]);
	    }
#endif
	}
      else
	{
	  /*
	    even though we skip this event as it could not make new gametes,
	    we need to account for the number of breakpoints that would have occured
	  */
	  //initialize a discrete distribution from 1 to nrecs that is a truncated Poisson
	  std::vector<double> breakpoint_dist;
	  for(unsigned i = 1 ; i <= NRECS ; ++i)
	    {
	      breakpoint_dist.push_back( gsl_ran_binomial_pdf(i,2.*littler,NRECS) );
	    }
	  lookup_t BPLOOKUP(gsl_ran_discrete_preproc (NRECS, &breakpoint_dist[0]));
	  NRECS_i = gsl_ran_discrete(r,BPLOOKUP.get())+1;
	  //gsl_ran_discrete_free(BPLOOKUP);
	}
      NRECS-=unsigned(NRECS_i);
#ifndef NDEBUG
      NRECS_DONE+=unsigned(NRECS_i);
#endif
    }
  assert(NRECS_DONE==doublecheck);
#ifndef NDEBUG
  unsigned sum=0;
  for(ibeg=gametes->begin();ibeg!=gametes->end();++ibeg)
    {
      sum+=ibeg->n;
    }
  assert(sum==twoN);
#endif
  //erase all gametes that have gone extinct...
  if(NEXTINCT)
    {
      vtype_iterator newend = std::remove_if(gametes->begin(),gametes->end(),
					     std::bind(n_is_zero(),std::placeholders::_1));
      gametes->erase(newend,gametes->end());
    }
  return NRECS;
}
}

#endif
