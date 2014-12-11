//  -*- C++ -*- 
#ifndef _MUTATION_TCC_
#define _MUTATION_TCC_

#include <type_traits>
#include <algorithm>
#include <numeric>

#include <gsl/gsl_randist.h>

#if defined(HAVE_BOOST_VECTOR) && !defined(USE_STANDARD_CONTAINERS)
#include <boost/container/vector.hpp>
#endif

#include <fwdpp/internal/mutation_internal.hpp>
#include <vector>

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
	      decltype(gametes->size()) ncurrent_classes = gametes->size();
	      decltype(gametes->begin()) ibeg;
	      unsigned NM=0;

	      for(decltype(ncurrent_classes) i=0;i<ncurrent_classes;++i)
		{
		  assert( i < std::numeric_limits<typename std::iterator_traits<decltype(ibeg)>::difference_type>::max() );
		  ibeg=(gametes->begin()+typename std::iterator_traits<decltype(ibeg)>::difference_type(i));
		  unsigned nmuts = gsl_ran_poisson(r,double(ibeg->n)*mu);
		  NM += nmuts;
		  
#if defined(HAVE_BOOST_VECTOR) && !defined(USE_STANDARD_CONTAINERS)
		  typedef boost::container::vector<unsigned> vu;
		  boost::container::vector<double> pm(ibeg->n,1./double(ibeg->n));
#else
		  typedef std::vector<unsigned> vu;
		  std::vector<double> pm(ibeg->n,1./double(ibeg->n));
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
		      fwdpp_internal::add_N_mutations(mmodel,mpolicy,*itr,mutations,new_gamete);
		      nmuts -= *itr;
		      gpolicy(new_gamete,gametes);
		      ibeg=(gametes->begin()+typename std::iterator_traits<decltype(ibeg)>::difference_type(i));
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
    if ( nm )
      {
	assert( g->n > 0 );
	g->n--;
	typename iterator_type::value_type ng( 1, g->mutations,g->smutations);
	fwdpp_internal::add_N_mutations(mmodel,mpolicy,nm,mutations,ng);
	return gpolicy(ng,gametes);
      }
    return g;
  }
}
#endif /* _MUTATION_TCC_ */
