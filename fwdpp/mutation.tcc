//  -*- C++ -*- 
#ifndef _MUTATION_TCC_
#define _MUTATION_TCC_

#include <type_traits>
#include <algorithm>
#include <numeric>

#include <gsl/gsl_randist.h>

#ifdef USE_STANDARD_CONTAINERS
#include <vector>
#else
#include <boost/container/vector.hpp>
#endif

namespace KTfwd
{
  template< typename gamete_type,
	    typename mutation_model,
	    typename gamete_insertion_policy,
	    typename mutation_insertion_policy,
	    typename vector_type_allocator,
	    typename list_type_allocator,
	    template<typename,typename> class vector_type,
	    template<typename,typename> class list_type>
  unsigned mutate(gsl_rng * r, 
		  vector_type<gamete_type,vector_type_allocator > * gametes, 
		  list_type<typename gamete_type::mutation_type,list_type_allocator > * mutations, 
		  const double & mu,
		  const mutation_model & mmodel,
		  const gamete_insertion_policy & gpolicy, 
		  const mutation_insertion_policy & mpolicy)
	    {
	      static_assert( std::is_base_of<mutation_base,typename gamete_type::mutation_type>::value,
                             "gamete_type::mutation_type must be derived from KTfwd::mutation_base" );
	      typedef gamete_base< typename gamete_type::mutation_type, 
				   typename gamete_type::mutation_list_type > gamete_base_type;
	      static_assert( std::is_base_of<gamete_base_type,gamete_type>::value ||
                             std::is_same<gamete_base_type,gamete_type>::value,
                             "gamete_type must be, or inherit from, KTfwd::gamete_base<mutation_type,mutation_list_type>" );
              static_assert( std::is_same<list_type<typename gamete_type::mutation_type,list_type_allocator >,
                                          typename gamete_type::mutation_list_type >::value,
                             "list_type<typename gamete_type::mutation_type,list_type_allocator > and gamete_type::mutation_list_type must be the same" );
	      unsigned ncurrent_classes = gametes->size();
	      typename vector_type<gamete_type,vector_type_allocator>::iterator ibeg;
	      typedef typename gamete_type::mutation_type mut_type;
	      typedef typename list_type<typename gamete_type::mutation_type,list_type_allocator >::iterator mut_itr;
	      unsigned NM=0;

	      for(unsigned i=0;i<ncurrent_classes;++i)
		{
		  ibeg=(gametes->begin()+i);
		  unsigned nmuts = gsl_ran_poisson(r,double(ibeg->n)*mu);
		  NM += nmuts;
		  
#ifdef USE_STANDARD_CONTAINERS
		  typedef std::vector<unsigned> vu;
		  std::vector<double> pm(ibeg->n,1./double(ibeg->n));
#else
		  typedef boost::container::vector<unsigned> vu;
		  boost::container::vector<double> pm(ibeg->n,1./double(ibeg->n));
#endif
		  vu nm(ibeg->n,0u);
		  gsl_ran_multinomial(r,ibeg->n,nmuts,&pm[0],&nm[0]);
		  
		  assert( std::accumulate(nm.begin(),nm.end(),0u) == nmuts );
		  
		  std::sort(nm.begin(),nm.end(),std::greater<unsigned>());
		  for(vu::const_iterator itr = nm.begin() ; 
		      itr < nm.end() && *itr>0 ; ++itr )
		    {
		      ibeg->n--;
		      gamete_type new_gamete( 1,ibeg->mutations,ibeg->smutations );
		      for(unsigned j=0;j<*itr;++j)
			{
			  nmuts--;
			  
			  //create a new mutant type to enter the population
			  mut_type new_mutant = mmodel(mutations);
			  
			  //insert the new mutant type into the list of mutations, and record the position of insertion
			  mut_itr mitr = mpolicy(new_mutant,mutations);
			  
			  if(mitr->neutral)
			    {
			      new_gamete.mutations.insert( (new_gamete.mutations.size()<=250) ? 
							   std::find_if(new_gamete.mutations.begin(),
									new_gamete.mutations.end(),
									[&](const typename gamete_type::mutation_list_type_iterator & mut_itr)
									{
									  return mut_itr->pos > mitr->pos;
									}) :
							   std::upper_bound(new_gamete.mutations.begin(),
									    new_gamete.mutations.end(),mitr->pos,
									    [](const double & __value,const mut_itr & __mut){ return __value < __mut->pos;}),
							   mitr );
			    }
			  else
			    {
			      new_gamete.smutations.insert( (new_gamete.smutations.size() <= 250) ? 
							    std::find_if(new_gamete.smutations.begin(),
									 new_gamete.smutations.end(),
									 [&](const typename gamete_type::mutation_list_type_iterator & mut_itr)
									 {
									   return mut_itr->pos > mitr->pos;
									 }) :
							    std::upper_bound(new_gamete.smutations.begin(),
									     new_gamete.smutations.end(),mitr->pos,
									     [](const double & __value,const mut_itr & __mut){ return __value < __mut->pos;}),
							    mitr );
			    }
			}
		      gpolicy(new_gamete,gametes);
		      ibeg=(gametes->begin()+i);
		    }
		  assert(nmuts==0);
		}
	      return NM;
	    }
  
  template< typename iterator_type,
	    typename mutation_model,
	    typename mutation_insertion_policy,
	    typename gamete_insertion_policy,
	    typename list_type_allocator,
	    typename list_type_allocator2,
	    template<typename,typename> class list_type,
	    template<typename,typename> class list_type2>
  iterator_type mutate_gamete( gsl_rng * r,
			       const double & mu, list_type< typename iterator_type::value_type,list_type_allocator > * gametes,
			       list_type2<typename iterator_type::value_type::mutation_type,list_type_allocator2 > * mutations, 
			       iterator_type &g,
			       const mutation_model &mmodel,
			       const mutation_insertion_policy & mpolicy,
			       const gamete_insertion_policy & gpolicy)
  {
    assert( g != gametes->end() );
    unsigned nm = gsl_ran_poisson(r,mu);
    typedef typename iterator_type::value_type::mutation_list_type_iterator mut_itr;
    if ( nm )
      {
	assert( g->n > 0 );
	g->n--;
	typename iterator_type::value_type ng( *g );
	ng.n = 1;
	for( unsigned i = 0 ; i < nm ; ++i )
	  {
	    mut_itr mitr = mpolicy(mmodel(mutations),mutations);
	    if( mitr->neutral )
	      {
		ng.mutations.insert((ng.mutations.size() <= 250) ? std::find_if(ng.mutations.begin(),
										ng.mutations.end(),
										[&](const mut_itr & __mut){return __mut->pos > mitr->pos;}) :
				    std::upper_bound(ng.mutations.begin(),
						     ng.mutations.end(),mitr->pos,
						     [](const double & __value,const mut_itr & __mut){ return __value < __mut->pos;}),
				    mitr);
	      }
	    else
	      {
		ng.smutations.insert((ng.smutations.size() <= 250) ? std::find_if(ng.smutations.begin(),
										  ng.smutations.end(),
										  [&](const mut_itr & __mut){return __mut->pos > mitr->pos;}) :
				     std::upper_bound(ng.smutations.begin(),
						      ng.smutations.end(),mitr->pos,
						      [](const double & __value,const mut_itr & __mut){ return __value < __mut->pos;}),
				     mitr);
	      }
	  }
	return gpolicy(ng,gametes);
      }
    return g;
  }
}
#endif /* _MUTATION_TCC_ */
