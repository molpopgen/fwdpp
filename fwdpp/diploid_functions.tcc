//  -*- C++ -*- 
#ifndef _DIPLOID_FUNCTIONS_TCC_
#define _DIPLOID_FUNCTIONS_TCC_

#include <utility>
#include <cassert>
#include <algorithm>
#include <iterator>
#include <functional>

#include <cmath>
#include <numeric>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/is_base_and_derived.hpp>
#include <boost/container/flat_map.hpp>
#include <boost/container/map.hpp>
#include <boost/container/vector.hpp>
#include <boost/pool/pool_alloc.hpp>
#include <fwdpp/util.hpp>
#include <fwdpp/insertion_policies.hpp>
#include <fwdpp/fwd_functional.hpp>

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
	      BOOST_STATIC_ASSERT( (boost::is_base_and_derived<mutation_base,typename gamete_type::mutation_type>::value) );
	      typedef gamete_base< typename gamete_type::mutation_type, typename gamete_type::mutation_list_type > gamete_base_type;
	      BOOST_STATIC_ASSERT( (boost::is_base_and_derived<gamete_base_type,
				    gamete_type>::value) || (boost::is_same<gamete_base_type,gamete_type >::value) );
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
  BOOST_STATIC_ASSERT( (boost::is_base_and_derived<mutation_base,typename gamete_type::mutation_type>::value) );
  typedef gamete_base< typename gamete_type::mutation_type, typename gamete_type::mutation_list_type > gamete_base_type;
  BOOST_STATIC_ASSERT( (boost::is_base_and_derived<gamete_base_type,
			gamete_type>::value) || (boost::is_same<gamete_base_type,gamete_type >::value) );
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
	      BOOST_STATIC_ASSERT( (boost::is_base_and_derived<mutation_base,typename gamete_type::mutation_type>::value) );
	      typedef gamete_base< typename gamete_type::mutation_type, typename gamete_type::mutation_list_type > gamete_base_type;
	      BOOST_STATIC_ASSERT( (boost::is_base_and_derived<gamete_base_type,
				    gamete_type>::value) || (boost::is_same<gamete_base_type,gamete_type >::value) );
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
	      BOOST_STATIC_ASSERT( (boost::is_base_and_derived<mutation_base,typename gamete_type::mutation_type>::value) );
	      typedef gamete_base< typename gamete_type::mutation_type, typename gamete_type::mutation_list_type > gamete_base_type;
	      BOOST_STATIC_ASSERT( (boost::is_base_and_derived<gamete_base_type,
				    gamete_type>::value) || (boost::is_same<gamete_base_type,gamete_type >::value) );
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
	      BOOST_STATIC_ASSERT( (boost::is_base_and_derived<mutation_base,typename gamete_type::mutation_type>::value) );
	      typedef gamete_base< typename gamete_type::mutation_type, typename gamete_type::mutation_list_type > gamete_base_type;
	      BOOST_STATIC_ASSERT( (boost::is_base_and_derived<gamete_base_type,
				    gamete_type>::value) || (boost::is_same<gamete_base_type,gamete_type >::value) );
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
	      BOOST_STATIC_ASSERT( (boost::is_base_and_derived<mutation_base,typename gamete_type::mutation_type>::value) );
	      typedef gamete_base< typename gamete_type::mutation_type, typename gamete_type::mutation_list_type > gamete_base_type;
	      BOOST_STATIC_ASSERT( (boost::is_base_and_derived<gamete_base_type,
				    gamete_type>::value) || (boost::is_same<gamete_base_type,gamete_type >::value) );
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
	      BOOST_STATIC_ASSERT( (boost::is_base_and_derived<mutation_base,typename gamete_type::mutation_type>::value) );
	      typedef gamete_base< typename gamete_type::mutation_type, typename gamete_type::mutation_list_type > gamete_base_type;
	      BOOST_STATIC_ASSERT( (boost::is_base_and_derived<gamete_base_type,
				    gamete_type>::value) || (boost::is_same<gamete_base_type,gamete_type >::value) );
	      
	      typedef typename gamete_type::mutation_container mcont;
	      typedef typename mcont::iterator mcont_iterator;
	      typedef typename vector_type<gamete_type, vector_type_allocator>::iterator vtype_iterator;
	      //"individual-based recombination"
	      //1. determine # recombinants in whole pop
	      //2. assign them randomly to pairs of gametes chosen from the gametes list
	      //3. if they are the same gamete, ignore, else do recombination thing.

	      //calc freq of homozygotes
  
	      double fAA=0.;
#ifndef USE_STANDARD_CONTAINERS
	      boost::container::vector<double> gcounts(gametes->size(),0);
#else
	      std::vector<double> gcounts(gametes->size(),0);
#endif
	      const double twoNsq = double(twoN)*double(twoN);
	      for(vtype_iterator i=gametes->begin() ;i<gametes->end();++i)
		{
		  //gcounts.push_back(i->n);
		  gcounts[i - gametes->begin()] = i->n;
		  fAA += ( (double(i->n)*double(i->n)*(1.-f))/twoNsq + (double(i->n)/twoN)*f );
		}
	      assert(gcounts.size()==gametes->size());
	      gsl_ran_discrete_t * lookup = gsl_ran_discrete_preproc(gcounts.size(),&gcounts[0]);
	      unsigned NRECS = unsigned(gsl_ran_poisson(r,(1.-fAA)*double(twoN)*littler));
#ifndef NDEBUG
	      unsigned doublecheck=NRECS;
	      unsigned NRECS_DONE = 0;
#endif
	      vtype_iterator ibeg,jbeg;
	      unsigned ncurrent_classes = gametes->size(),NEXTINCT=0,NRECS_i;
	      size_t ith,jth;
	      while(NRECS > 0)
		{
		  ith = gsl_ran_discrete(r,lookup);
		  while(gcounts[ith]==0)
		    {
		      ith = gsl_ran_discrete(r,lookup);
		    }
		  jth = gsl_ran_discrete(r,lookup);
		  while(ith==jth || gcounts[jth]==0)
		    {
		      jth = gsl_ran_discrete(r,lookup);
		    }
		  ibeg = (gametes->begin()+ith);
		  jbeg = (gametes->begin()+jth);
		  assert( ibeg != gametes->end() &&
			  jbeg != gametes->end() );
		  assert(ibeg->n >= 1);
		  assert(jbeg->n >= 1);
		  unsigned nm1=ibeg->mutations.size()+ibeg->smutations.size(),
		    nm2=jbeg->mutations.size()+jbeg->smutations.size();
		  //if one gametes carries 0 mutations, and the other carries 1,
		  //it is impossible for recombination to generate any new types,
		  //so we skip those cases...
		  if(!(std::min(nm1,nm2)==0 && std::max(nm1,nm2)==1))
		    {
		      //initialize a discrete distribution from 1 to nrecs that is a truncated Poisson
#ifndef USE_STANDARD_CONTAINERS
		      boost::container::vector<double> breakpoint_dist;
#else
		      std::vector<double> breakpoint_dist;
#endif
		      for(unsigned i = 1 ; i <= NRECS ; ++i)
		      	{
		      	  breakpoint_dist.push_back( gsl_ran_poisson_pdf(i,2.*littler) );
		      	}
		      gsl_ran_discrete_t * BPLOOKUP = gsl_ran_discrete_preproc (NRECS, &breakpoint_dist[0]);
		      NRECS_i = gsl_ran_discrete(r,BPLOOKUP)+1;
		      gsl_ran_discrete_free(BPLOOKUP);
		      //assign breakpoints
#ifndef USE_STANDARD_CONTAINERS
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
		      NEXTINCT += (ibeg->n==0) ? 1 : 0;
		      NEXTINCT += (jbeg->n==0) ? 1 : 0;
		      gamete_type new_gamete1(1, mcont(),mcont()),new_gamete2(new_gamete1);
		      new_gamete1.mutations.reserve(ibeg->mutations.size()+jbeg->mutations.size());
		      new_gamete1.smutations.reserve(ibeg->smutations.size()+jbeg->smutations.size());
		      new_gamete2.mutations.reserve(ibeg->mutations.size()+jbeg->mutations.size());
		      new_gamete2.smutations.reserve(ibeg->smutations.size()+jbeg->smutations.size());
		      short SWITCH_I = 0,SWITCH_J = 0, SWITCH_Is = 0, SWITCH_Js = 0;
		      size_t dummy = 0;
		      mcont_iterator itr = ibeg->mutations.begin(),
			jtr = jbeg->mutations.begin(),
			itr_s = ibeg->smutations.begin(),
			jtr_s = jbeg->smutations.begin();
		      //pointer arithmetic over a range of pointers.  apologies...
		      //typename gamete_type::mutation_container::iterator itr2;
		      for( ; dummy < pos.size(); ++dummy)
			{
			  //iterate over neutral mutations from parent i
			  for( ; itr < ibeg->mutations.end() && (*itr)->pos < pos[dummy] ;++itr)
			    {
			      switch(SWITCH_I)
				{
				case 0:
#ifndef NDEBUG
				  for(unsigned t = 0 ; t < new_gamete1.mutations.size() ; ++t )
				    {
				      assert(! mutation_at_pos()(*new_gamete1.mutations[t],(*itr)->pos) );
				    }
#endif
				  new_gamete1.mutations.push_back(*itr);
				  break;
				case 1:
#ifndef NDEBUG
				  for(unsigned t = 0 ; t < new_gamete2.mutations.size() ; ++t )
				    {
				      assert(! mutation_at_pos()(*new_gamete2.mutations[t],(*itr)->pos) );
				    }
#endif
				  new_gamete2.mutations.push_back(*itr);
				  break;
				}
			    }
			  SWITCH_I=!SWITCH_I;

			  //iterate over selected mutations from parent i
			  for( ; itr_s < ibeg->smutations.end() && (*itr_s)->pos < pos[dummy] ;++itr_s)
			    {
			      switch(SWITCH_Is)
				{
				case 0:
#ifndef NDEBUG
				  for(unsigned t = 0 ; t < new_gamete1.smutations.size() ; ++t )
				    {
				      assert(! mutation_at_pos()(*new_gamete1.smutations[t],(*itr_s)->pos) );
				    }
#endif
				  new_gamete1.smutations.push_back(*itr_s);
				  break;
				case 1:
#ifndef NDEBUG
				  for(unsigned t = 0 ; t < new_gamete2.smutations.size() ; ++t )
				    {
				      assert(! mutation_at_pos()(*new_gamete2.smutations[t],(*itr_s)->pos) );
				    }
#endif
				  new_gamete2.smutations.push_back(*itr_s);
				  break;
				}
			    }
			  SWITCH_Is=!SWITCH_Is;

			  //iterate over neutral mutations from parent j
			  for( ; jtr < jbeg->mutations.end() && (*jtr)->pos < pos[dummy] ;++jtr)
			    {
			      switch(SWITCH_J)
				{
				case 1:
#ifndef NDEBUG
				  for(unsigned t = 0 ; t < new_gamete1.mutations.size() ; ++t )
				    {
				      assert(! mutation_at_pos()(*new_gamete1.mutations[t],(*jtr)->pos) );
				    }
#endif
				  new_gamete1.mutations.push_back(*jtr);
				  break;
				case 0:
#ifndef NDEBUG
				  for(unsigned t = 0 ; t < new_gamete2.mutations.size() ; ++t )
				    {
				      assert(! mutation_at_pos()(*new_gamete2.mutations[t],(*jtr)->pos) );
				    }
#endif
				  new_gamete2.mutations.push_back(*jtr);
				  break;
				}
			    }
			  SWITCH_J=!SWITCH_J;

 			  //iterate over selected mutations from parent j
			  for( ; jtr_s < jbeg->smutations.end() && (*jtr_s)->pos < pos[dummy] ;++jtr_s)
			    {
			      switch(SWITCH_Js)
				{
				case 1:
#ifndef NDEBUG
				  for(unsigned t = 0 ; t < new_gamete1.smutations.size() ; ++t )
				    {
				      assert(! mutation_at_pos()(*new_gamete1.smutations[t],(*jtr_s)->pos) );
				    }
#endif
				  new_gamete1.smutations.push_back(*jtr_s);
				  break;
				case 0:
#ifndef NDEBUG
				  for(unsigned t = 0 ; t < new_gamete2.smutations.size() ; ++t )
				    {
				      assert(! mutation_at_pos()(*new_gamete2.smutations[t],(*jtr_s)->pos) );
				    }
#endif
				  new_gamete2.smutations.push_back(*jtr_s);
				  break;
				}
			    }
			  SWITCH_Js=!SWITCH_Js;
			}
#ifndef NDEBUG
		      unsigned __nm1 = new_gamete1.mutations.size()+new_gamete1.smutations.size(),
			__nm2 = new_gamete2.mutations.size()+new_gamete2.smutations.size();
		      assert(__nm1+__nm2 == nm1+nm2);
#endif

		      std::sort(new_gamete1.mutations.begin(),new_gamete1.mutations.end(),fake_less());
		      std::sort(new_gamete1.smutations.begin(),new_gamete1.smutations.end(),fake_less());
		      std::sort(new_gamete2.mutations.begin(),new_gamete2.mutations.end(),fake_less());
		      std::sort(new_gamete2.smutations.begin(),new_gamete2.smutations.end(),fake_less());

		      vtype_iterator newpos = update_if_exists_insert(new_gamete1,gametes);
		      assert(newpos->n <= twoN);
		      size_t dist = newpos - gametes->begin();
		      if( dist < gcounts.size() )
			{
			  gcounts[dist]++;
			  assert( gcounts[dist] <= twoN );
			}
		      newpos = update_if_exists_insert(new_gamete2,gametes);
		      assert(newpos->n <= twoN);
		      dist = newpos - gametes->begin();
		      if( dist < gcounts.size() )
			{
			  gcounts[dist]++;
			  assert( gcounts[dist] <= twoN );
			}
#ifndef NDEBUG
		      dummy=0;
		      for( vtype_iterator test = gametes->begin();test!=gametes->begin()+ncurrent_classes ;++test,++dummy )
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
		      gsl_ran_discrete_t * BPLOOKUP = gsl_ran_discrete_preproc (NRECS, &breakpoint_dist[0]);
		      NRECS_i = gsl_ran_discrete(r,BPLOOKUP)+1;
		      gsl_ran_discrete_free(BPLOOKUP);
		    }
		  NRECS-=NRECS_i;
#ifndef NDEBUG
		  NRECS_DONE+=NRECS_i;
#endif
		}
	      gsl_ran_discrete_free(lookup);
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
							 boost::bind(n_is_zero(),_1));
		  gametes->erase(newend,gametes->end());
		}
	      return NRECS;
	    }

  //CODE FOR INDIVIDUAL-BASED FORWARD SIMULATIONS IS BELOW

  //single deme, constant N
  template< typename gamete_type,
	    typename gamete_list_type_allocator,
	    typename mutation_list_type_allocator,
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
		 diploid_vector_type<std::pair<typename gamete_list_type< gamete_type,gamete_list_type_allocator >::iterator,
					       typename gamete_list_type< gamete_type,gamete_list_type_allocator >::iterator>,
				     diploid_vector_type_allocator> * diploids,
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
		 diploid_vector_type<std::pair<typename gamete_list_type< gamete_type,gamete_list_type_allocator >::iterator,
					       typename gamete_list_type< gamete_type,gamete_list_type_allocator >::iterator>,
				     diploid_vector_type_allocator> * diploids,
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
    
    for( typename mutation_list_type<typename gamete_type::mutation_type,mutation_list_type_allocator>::iterator itr = mutations->begin() ; 
	 itr != mutations->end() ; ++itr )
      {
	itr->n = 0;
      }
	      
    typedef diploid_vector_type<std::pair<typename gamete_list_type< gamete_type,gamete_list_type_allocator >::iterator,
					  typename gamete_list_type< gamete_type,gamete_list_type_allocator >::iterator>,
				diploid_vector_type_allocator> dipctr;
    
    std::vector<double> fitnesses(diploids->size());
    double wbar = 0.;
    
    typename dipctr::iterator dptr = diploids->begin();
    for( unsigned i = 0 ; i < N_curr ; ++i )
      {
	(dptr+i)->first->n = 0;
	(dptr+i)->second->n = 0;
	fitnesses[i] = ff((dptr+i)->first,(dptr+i)->second);
	wbar += fitnesses[i];
      }
    wbar /= double(diploids->size());
    
#ifndef NDEBUG
    for(typename gamete_list_type<gamete_type,gamete_list_type_allocator >::iterator itr = gametes->begin() ; itr != gametes->end() ;++itr)
      {
	assert(itr->n==0);
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
	
	p1g1 = (pptr+p1)->first;
	p1g2 = (pptr+p1)->second;
	p2g1 = (pptr+p2)->first;
	p2g2 = (pptr+p2)->second;
	
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
	
	adjust_mutation_counts((dptr+i)->first,1);
	adjust_mutation_counts((dptr+i)->second,1);
	
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
    gametes->remove_if(boost::bind(n_is_zero(),_1));
    for( typename gamete_list_type<gamete_type,gamete_list_type_allocator >::iterator itr = gametes->begin() ; 
	 itr != gametes->end() ; ++itr )
      {
	itr->mutations.erase( std::remove_if(itr->mutations.begin(),itr->mutations.end(),mp),itr->mutations.end() );
	itr->smutations.erase( std::remove_if(itr->smutations.begin(),itr->smutations.end(),mp),itr->smutations.end() );
      }
    assert(check_sum(gametes,2*N_next));
    gsl_ran_discrete_free(lookup);
    return wbar;
  }
  
  //Metapopulation version of sample_diploid for individual-based simulations and constant N
  template< typename gamete_type,
	    typename metapop_gamete_vector_type_allocator,
	    typename metapop_diploid_vector_type_allocator,
	    typename gamete_list_type_allocator,
	    typename mutation_list_type_allocator,
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
	    template<typename,typename> class metapop_gamete_vector_type,
	    template<typename,typename> class metapop_diploid_vector_type>
  std::vector< double >
  sample_diploid(gsl_rng * r,
		 metapop_gamete_vector_type < gamete_list_type<gamete_type,gamete_list_type_allocator > ,
		 metapop_gamete_vector_type_allocator > * metapop,
		 metapop_diploid_vector_type < diploid_vector_type<std::pair<typename gamete_list_type< gamete_type,gamete_list_type_allocator >::iterator,
									     typename gamete_list_type< gamete_type,gamete_list_type_allocator >::iterator>,
								   diploid_vector_type_allocator>,
					       metapop_diploid_vector_type_allocator > * diploids,
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
	    typename metapop_gamete_vector_type_allocator,
	    typename metapop_diploid_vector_type_allocator,
	    typename gamete_list_type_allocator,
	    typename mutation_list_type_allocator,
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
	    template<typename,typename> class metapop_gamete_vector_type,
	    template<typename,typename> class metapop_diploid_vector_type>
  std::vector< double >
  sample_diploid(gsl_rng * r,
		 metapop_gamete_vector_type < gamete_list_type<gamete_type,gamete_list_type_allocator > ,
		 metapop_gamete_vector_type_allocator > * metapop,
		 metapop_diploid_vector_type < diploid_vector_type<std::pair<typename gamete_list_type< gamete_type,gamete_list_type_allocator >::iterator,
									     typename gamete_list_type< gamete_type,gamete_list_type_allocator >::iterator>,
								   diploid_vector_type_allocator>,
					       metapop_diploid_vector_type_allocator > * diploids,
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
	      assert( metapop->size() == diploids->size() );
	      
	      for( typename mutation_list_type<typename gamete_type::mutation_type,mutation_list_type_allocator >::iterator mptr = mutations->begin() ; 
		   mptr != mutations->end() ; ++mptr )
		{
		  mptr->n = 0;
		}
	      
	      typedef metapop_gamete_vector_type < gamete_list_type<gamete_type,gamete_list_type_allocator > ,
						   metapop_gamete_vector_type_allocator > pop_ctr;

	      typedef gamete_list_type<gamete_type,gamete_list_type_allocator > gamete_ctr;
	      
	      //we need tpop b/c saying &*pop_ptr below will result in passing a const pointer on at least some compilers (e.g., like mine, which seems lame)
	      gamete_ctr * tpop;
	      
	      typedef metapop_diploid_vector_type < diploid_vector_type<std::pair<typename gamete_ctr::iterator,
										  typename gamete_ctr::iterator>,
									diploid_vector_type_allocator>,
						    metapop_diploid_vector_type_allocator > diploid_ctr;
	      
	      typedef diploid_vector_type<std::pair<typename gamete_ctr::iterator,
						    typename gamete_ctr::iterator>,
					  diploid_vector_type_allocator> demedips;
  
	      //get the fitnesses for each diploid in each deme and make the lookup table of parental fitnesses
	      std::vector<gsl_ran_discrete_t *> lookups(metapop->size());
	      std::vector<double> wbars(metapop->size(),0);
	      size_t popindex = 0;

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
	      
	      for( typename diploid_ctr::iterator dptr = diploids->begin() ; dptr != diploids->end() ; ++dptr, ++popindex )
		{
		  unsigned demesize = *(N_curr+popindex);
		  assert( demesize == dptr->size() );
		  size_t ith_dip = 0;
		  for( typename demedips::iterator gptr = dptr->begin() ; 
		       gptr != dptr->end() ; ++gptr,++ith_dip )
		    {
		      fitnesses[ith_dip] = ffs[popindex](gptr->first,gptr->second);
		      wbars[popindex]+=fitnesses[ith_dip];
		      gptr->first->n = 0;
		      gptr->second->n = 0;
		    }
		  wbars[popindex] /= double( demesize );
		  lookups[popindex]=gsl_ran_discrete_preproc(demesize,fitnesses);
		}
	      delete [] fitnesses;
	      
	      assert(lookups.size() == diploids->size());
	      //copy diploids into temporary parents
	      diploid_ctr parents(*diploids);
	      
	      //Update the diploids
	      popindex = 0;
	      unsigned NREC=0;
	      typename pop_ctr::iterator pop_ptr = metapop->begin();
	      typename gamete_ctr::iterator p1g1,p1g2,p2g1,p2g2;
	      for( typename diploid_ctr::iterator ptr = diploids->begin() ; ptr != diploids->end() ; ++ptr, ++pop_ptr,++popindex )
		{
		  unsigned demesize =*(N_next+popindex);
		  if( demesize != *(N_curr+popindex) )
		    {
		      ptr->resize(demesize);
		    }
		  const typename demedips::iterator dptr = ptr->begin();
		  //typename demedips::iterator pptr=(parents.begin()+popindex)->begin(),pptr2;
		  for( unsigned i = 0 ; i < demesize ; ++i )
		    {
		      //figure out if parent 1 is migrant or not
		      size_t deme_first_parent = mig(popindex),deme_other_parent=popindex;
		      typename demedips::iterator pptr=(parents.begin()+deme_first_parent)->begin();
		      size_t p1 = gsl_ran_discrete(r,lookups[deme_first_parent]),p2;
		      
		      tpop = &*pop_ptr;
		      
		      if( popindex == deme_first_parent )
			//not migrant
			{
			  p1g1 = (pptr+p1)->first;
			  p1g2 = (pptr+p1)->second;
			}
		      else
			//migrant                                                                                                    
			{
			  p1g1 = insert_if_not_found( *((pptr+p1)->first),tpop );
			  p1g2 = insert_if_not_found( *((pptr+p1)->second),tpop );
			}
		      
		      NREC += rec_pol(p1g1,p1g2,tpop);
		      (dptr+i)->first = (gsl_rng_uniform(r) <= 0.5) ? p1g1 : p1g2;
		      
		      /*
			If the individual is not inbred, then we pick a 
			deme from the migration policy.
			
			A migration policy takes the current deme (popindex) as
			an argument.  It returns popindex if there is no migration,
			else it returns the index of the deme of a migrant parent
		      */
		      typename demedips::iterator pptr2=(parents.begin()+deme_other_parent)->end();
		      if( f != NULL && gsl_rng_uniform(r) <= *(f + popindex ) ) //individual is inbred
			{
			  pptr2=(parents.begin()+popindex)->begin();
			  p2=p1;
			}
		      else
			{
			  deme_other_parent = mig(popindex);
			  assert( deme_other_parent < diploids->size() );
			  pptr2 = (parents.begin() + deme_other_parent)->begin();
			  p2 = gsl_ran_discrete(r,lookups[deme_other_parent]);
			  assert( (pptr2+p2) < (parents.begin() + deme_other_parent)->end() );
			}
		      assert( pptr2 != (parents.begin() + deme_other_parent)->end() );
		      if(deme_other_parent == popindex)
			{
			  p2g1 = (pptr2+p2)->first;
			  p2g2 = (pptr2+p2)->second;
			}
		      else
			{
			  //We may need to put p2's gametes into the pop pointed to by pop_ptr
			  p2g1 = insert_if_not_found( *((pptr2+p2)->first),tpop);
			  p2g2 = insert_if_not_found( *((pptr2+p2)->second),tpop );
			}
		      
		      NREC += rec_pol( p2g1,p2g2, tpop );
		      (dptr+i)->second = (gsl_rng_uniform(r) <= 0.5) ? p2g1 : p2g2;
		      assert( std::find( (pop_ptr)->begin(), (pop_ptr)->end(), *( (dptr+i)->second ) )
			      != (pop_ptr)->end() );
		      
		      (dptr+i)->first->n++;
		      assert((dptr+i)->first->n <= 2*demesize);
		      (dptr+i)->second->n++;
		      assert((dptr+i)->second->n <= 2*demesize);
		      
		      adjust_mutation_counts((dptr+i)->first,1);
		      adjust_mutation_counts((dptr+i)->second,1);
		      
		      //now, add new mutations
		      (dptr+i)->first = mutate_gamete(r,mu,&*(pop_ptr),mutations,(dptr+i)->first,mmodel,mpolicy,gpolicy_mut);
		      (dptr+i)->second = mutate_gamete(r,mu,&*(pop_ptr),mutations,(dptr+i)->second,mmodel,mpolicy,gpolicy_mut);
		    }
		}
	      
	      //get rid of extinct stuff, etc.
	      for(typename pop_ctr::iterator ptr = metapop->begin() ; ptr != metapop->end() ; ++ptr)
		{
		  ptr->remove_if(boost::bind(n_is_zero(),_1));
		  for (typename gamete_ctr::iterator gptr = ptr->begin() ; gptr != ptr->end() ; ++gptr )
		    {
		      gptr->mutations.erase( std::remove_if(gptr->mutations.begin(),gptr->mutations.end(),mp), gptr->mutations.end() );
		      gptr->smutations.erase( std::remove_if(gptr->smutations.begin(),gptr->smutations.end(),mp), gptr->smutations.end() );
		    }
		}
	      
	      for(unsigned i = 0 ; i < lookups.size() ; ++i )
		{
		  gsl_ran_discrete_free( lookups[i] );
		}
	      return wbars;
	    }
  
  //recombination for individual-based simulation
  template< typename iterator_type,
	    typename gamete_insertion_policy,
	    typename recombination_map,
	    typename list_type_allocator,
	    template<typename,typename> class list_type>
  unsigned recombine_gametes( gsl_rng * r,
			      const double & littler,
			      list_type< typename iterator_type::value_type,list_type_allocator > * gametes,
			      iterator_type & g1,
			      iterator_type & g2,
			      const recombination_map & mf,
			      const gamete_insertion_policy & gpolicy)
  {
    assert( g1 != gametes->end() );
    assert( g2 != gametes->end() );
    
    //Identify cases where recombination cannot result in changed gametes, and get out quick
    if(g1 == g2 ) return 0;
    unsigned nm1=g1->mutations.size()+g1->smutations.size();
    unsigned nm2=g2->mutations.size()+g2->smutations.size();
    if((std::min(nm1,nm2)==0 && std::max(nm1,nm2)==1)) return 0;
    
    unsigned nbreaks = (littler > 0) ? gsl_ran_poisson(r,littler) : 0.;
    
    if( nbreaks )
      {
#ifndef USE_STANDARD_CONTAINERS
	boost::container::vector<double> pos;
#else
	std::vector<double> pos;
#endif
	for(unsigned i = 0 ; i < nbreaks ; ++i)
	  {
	    pos.push_back(mf());
	  }
	std::sort(pos.begin(),pos.end());
	pos.push_back(std::numeric_limits<double>::max());
	
	typename iterator_type::value_type new_gamete1(0u,
						       typename iterator_type::value_type::mutation_container(),
						       typename iterator_type::value_type::mutation_container()),
	  new_gamete2(new_gamete1);
	new_gamete1.mutations.reserve(g1->mutations.size()+g2->mutations.size());
	new_gamete1.smutations.reserve(g1->smutations.size()+g2->smutations.size());
	new_gamete2.mutations.reserve(g1->mutations.size()+g2->mutations.size());
	new_gamete2.smutations.reserve(g1->smutations.size()+g2->smutations.size());
	short SWITCH_I = 0,SWITCH_J = 0, SWITCH_Is = 0, SWITCH_Js = 0;
	size_t dummy = 0;
	typename iterator_type::value_type::mcont_iterator itr = g1->mutations.begin(),
	  jtr = g2->mutations.begin(),
	  itr_s = g1->smutations.begin(),
	  jtr_s = g2->smutations.begin();
	//pointer arithmetic over a range of pointers.  apologies...
	//typename gamete_type::mutation_container::iterator itr2;
	for( ; dummy < pos.size(); ++dummy)
	  {
	    //iterate over neutral mutations from parent i
	    for( ; itr < g1->mutations.end() && (*itr)->pos < pos[dummy] ;++itr)
	      {
		switch(SWITCH_I)
		  {
		  case 0:
#ifndef NDEBUG
		    for(unsigned t = 0 ; t < new_gamete1.mutations.size() ; ++t )
		      {
			assert(! mutation_at_pos()(*new_gamete1.mutations[t],(*itr)->pos) );
		      }
#endif
		    new_gamete1.mutations.push_back(*itr);
		    break;
		  case 1:
#ifndef NDEBUG
		    for(unsigned t = 0 ; t < new_gamete2.mutations.size() ; ++t )
		      {
			assert(! mutation_at_pos()(*new_gamete2.mutations[t],(*itr)->pos) );
		      }
#endif
		    new_gamete2.mutations.push_back(*itr);
		    break;
		  }
	      }
	    SWITCH_I=!SWITCH_I;
	    
	    //iterate over selected mutations from parent i
	    for( ; itr_s < g1->smutations.end() && (*itr_s)->pos < pos[dummy] ;++itr_s)
	      {
		switch(SWITCH_Is)
		  {
		  case 0:
#ifndef NDEBUG
		    for(unsigned t = 0 ; t < new_gamete1.smutations.size() ; ++t )
		      {
			assert(! mutation_at_pos()(*new_gamete1.smutations[t],(*itr_s)->pos) );
		      }
#endif
		    new_gamete1.smutations.push_back(*itr_s);
		    break;
		  case 1:
#ifndef NDEBUG
		    for(unsigned t = 0 ; t < new_gamete2.smutations.size() ; ++t )
		      {
			assert(! mutation_at_pos()(*new_gamete2.smutations[t],(*itr_s)->pos) );
		      }
#endif
		    new_gamete2.smutations.push_back(*itr_s);
		    break;
		  }
	      }
	    SWITCH_Is=!SWITCH_Is;
	    
	    //iterate over neutral mutations from parent j
	    for( ; jtr < g2->mutations.end() && (*jtr)->pos < pos[dummy] ;++jtr)
	      {
		switch(SWITCH_J)
		  {
		  case 1:
#ifndef NDEBUG
		  for(unsigned t = 0 ; t < new_gamete1.mutations.size() ; ++t )
		    {
		      assert(! mutation_at_pos()(*new_gamete1.mutations[t],(*jtr)->pos) );
		    }
#endif
		  new_gamete1.mutations.push_back(*jtr);
		  break;
		case 0:
#ifndef NDEBUG
		  for(unsigned t = 0 ; t < new_gamete2.mutations.size() ; ++t )
		    {
		      assert(! mutation_at_pos()(*new_gamete2.mutations[t],(*jtr)->pos) );
		    }
#endif
		  new_gamete2.mutations.push_back(*jtr);
		  break;
		}
	    }
	  SWITCH_J=!SWITCH_J;
	  
	  //iterate over selected mutations from parent j
	  for( ; jtr_s < g2->smutations.end() && (*jtr_s)->pos < pos[dummy] ;++jtr_s)
	    {
	      switch(SWITCH_Js)
		{
		case 1:
#ifndef NDEBUG
		  for(unsigned t = 0 ; t < new_gamete1.smutations.size() ; ++t )
		    {
		      assert(! mutation_at_pos()(*new_gamete1.smutations[t],(*jtr_s)->pos) );
		    }
#endif
		  new_gamete1.smutations.push_back(*jtr_s);
		  break;
		case 0:
#ifndef NDEBUG
		  for(unsigned t = 0 ; t < new_gamete2.smutations.size() ; ++t )
		    {
		      assert(! mutation_at_pos()(*new_gamete2.smutations[t],(*jtr_s)->pos) );
		    }
#endif
		  new_gamete2.smutations.push_back(*jtr_s);
		  break;
		}
	    }
	  SWITCH_Js=!SWITCH_Js;
	}
	std::sort(new_gamete1.mutations.begin(),new_gamete1.mutations.end(),fake_less());
	std::sort(new_gamete1.smutations.begin(),new_gamete1.smutations.end(),fake_less());
	std::sort(new_gamete2.mutations.begin(),new_gamete2.mutations.end(),fake_less());
	std::sort(new_gamete2.smutations.begin(),new_gamete2.smutations.end(),fake_less());
	g1 = gpolicy(new_gamete1,gametes);
	g2 = gpolicy(new_gamete2,gametes);
      }
    return nbreaks;
  }
  
}//ns KTfwd
#endif /* _DIPLOID_FUNCTIONS_TCC_ */
