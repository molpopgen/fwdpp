/// \file util.hpp
#ifndef _UTIL_HPP_
#define _UTIL_HPP_

#include <fwdpp/forward_types.hpp>
#include <fwdpp/fwd_functional.hpp>
#include <fwdpp/internal/gsl_discrete.hpp>
#include <set>
#include <map>
#include <type_traits>
#include <algorithm>
#include <functional>
#include <cassert>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#if defined(HAVE_BOOST_VECTOR) && !defined(USE_STANDARD_CONTAINERS)
#include <boost/container/vector.hpp>
#else
#include <vector>
#endif


namespace KTfwd
{
  template< typename gamete_type,
	    typename gamete_allocator_type,
	    typename gamete_allocator_type2,
	    typename mutation_type,
	    typename mutation_allocator_type,
	    typename mutation_allocator_type2,
	    template<typename,typename> class source_gamete_container,
	    template<typename,typename> class source_mut_container,
	    template<typename,typename> class dest_gamete_container,
	    template<typename,typename> class dest_mut_container>
  void valid_copy( const source_gamete_container<gamete_type,gamete_allocator_type> & gametes,
		   const source_mut_container<mutation_type,mutation_allocator_type> & mutations,
		   dest_gamete_container<gamete_type,gamete_allocator_type2> & gametes_destination,
		   dest_mut_container<mutation_type,mutation_allocator_type2> & mutations_destination )
  /*!
    If you ever need to store (and later restore) the state of the population, a naive copy operation
    is not sufficient, because of all the pointers from the gametes container to elements
    of the mutations container.  Use this function instead.

    \note Only works for the case of unique mutation positions!!!
  */
  {
    static_assert( std::is_same<typename gamete_type::mutation_type,mutation_type>::value,
                   "gamete_type::mutation type and mutation_type must be the same" );
    typedef typename source_gamete_container<gamete_type,gamete_allocator_type>::const_iterator giterator;
    typedef typename source_mut_container<mutation_type,mutation_allocator_type>::iterator literator;  
    typedef typename source_mut_container<mutation_type,mutation_allocator_type>::const_iterator cliterator;
    typedef typename gamete_type::mutation_container::const_iterator gciterator;
    gametes_destination.clear();
    mutations_destination.clear();
    //copying the mutations is trivial
    std::map<double,literator> mutlookup;
    for( cliterator i = mutations.begin();
	 i!=mutations.end();++i)
      {
	literator j = mutations_destination.insert(mutations_destination.end(),*i);
	mutlookup[j->pos]=j;
      }
    
    for( giterator i = gametes.begin() ; i != gametes.end() ; ++i )
      {
	//copy construct so that all public, etc., data
	//are properly initialized
	gamete_type new_gamete(*i);//gametes[i]);
	new_gamete.mutations.clear();
	new_gamete.smutations.clear();
	for(gciterator itr = i->mutations.begin() ; 
	    itr != i->mutations.end() ; ++itr)
	  {
	    new_gamete.mutations.push_back( mutlookup[(*itr)->pos] );
	  }
	for(gciterator itr = i->smutations.begin() ; 
	    itr != i->smutations.end() ; ++itr)
	  {
	    new_gamete.smutations.push_back( mutlookup[(*itr)->pos] );
	  }
	gametes_destination.insert(gametes_destination.end(),new_gamete);
      }
  }

  /* \brief Sets mutation::checked to false
   */
  template<typename mutation_type,
	   typename list_type_allocator,
	   template <typename,typename> class list_type>
  void uncheck( list_type<mutation_type,list_type_allocator> * mutations )
  {
    static_assert( std::is_base_of<mutation_base,mutation_type>::value,
                   "mutation_type must be derived from KTfwd::mutation_base" );
    std::for_each(mutations->begin(),mutations->end(),[](mutation_type & __m){__m.checked=false;});
  }
  
  /*! \brief Remove mutations from population
    Removes mutations that are lost.
    \example diploid.cc
  */
  template<typename mutation_type,
	   typename list_type_allocator,
	   template <typename,typename> class list_type >
  void remove_lost( list_type<mutation_type,list_type_allocator> * mutations )
  {
    static_assert( std::is_base_of<mutation_base,mutation_type>::value,
                   "mutation_type must be derived from KTfwd::mutation_base" );
    typename list_type<mutation_type,list_type_allocator>::iterator i = mutations->begin(),
      temp;
    
    while(i != mutations->end())
      {
	i->checked = false;
	if( i->n == 0 )
	  {
	    temp = i;
	    ++i;
	    mutations->erase(temp);
	    //i=mutations->begin();
	  }
	else
	  {
	    ++i;
	  }
      }
  }

 /*! \brief Remove mutations from population
    Removes mutations that are lost.
    \note: lookup must be compatible with lookup->erase(lookup->find(double))
    \example diploid.cc
   */
  template<typename mutation_type,
	   typename list_type_allocator,
	   template <typename,typename> class list_type,
	   typename mutation_lookup_table>
  void remove_lost( list_type<mutation_type,list_type_allocator> * mutations, mutation_lookup_table * lookup )
  {
    static_assert( std::is_base_of<mutation_base,mutation_type>::value,
                   "mutation_type must be derived from KTfwd::mutation_base" );
    typename list_type<mutation_type,list_type_allocator>::iterator i = mutations->begin(),
      temp;
    
    while(i != mutations->end())
      {
	i->checked = false;
	if( i->n == 0 )
	  {
	    lookup->erase(lookup->find(i->pos));
	    temp=i;
	    ++i;
	    mutations->erase(temp);
	  }
	else
	  {
	    ++i;
	  }
      }
  }

  /*! \brief Remove mutations from population
    Removes mutations that are fixed or lost.
    \example diploid.cc
   */
  template<typename mutation_type,
	   typename vector_type_allocator1,
	   typename vector_type_allocator2,
	   typename list_type_allocator,
	   template <typename,typename> class vector_type,
	   template <typename,typename> class list_type >
  void remove_fixed_lost( list_type<mutation_type,list_type_allocator> * mutations, 
			  vector_type<mutation_type,vector_type_allocator1> * fixations, 
			  vector_type<unsigned,vector_type_allocator2> * fixation_times,
			  const unsigned & generation,const unsigned & twoN)
  {
    static_assert( std::is_base_of<mutation_base,mutation_type>::value,
		   "mutation_type must be derived from KTfwd::mutation_base" );
    typename list_type<mutation_type,list_type_allocator>::iterator i = mutations->begin(),temp;
    
    while(i != mutations->end())
      {
	assert(i->n <= twoN);			
	i->checked = false;
	if(i->n==twoN )
	  {
	    fixations->push_back(*i);
	    fixation_times->push_back(generation);
	  }
	if( i->n == 0 || i->n == twoN )
	  {
	    temp=i;
	    ++i;
	    mutations->erase(temp);
	  }
	else
	  {
	    ++i;
	  }
      }
  }

  /*! \brief Remove mutations from population
    Removes mutations that are fixed or lost.
    \note: lookup must be compatible with lookup->erase(lookup->find(double))
    \example diploid.cc
   */
  template<typename mutation_type,
	   typename vector_type_allocator1,
	   typename vector_type_allocator2,
	   typename list_type_allocator,
	   template <typename,typename> class vector_type,
	   template <typename,typename> class list_type,
	   typename mutation_lookup_table>
  void remove_fixed_lost( list_type<mutation_type,list_type_allocator> * mutations, 
			  vector_type<mutation_type,vector_type_allocator1> * fixations, 
			  vector_type<unsigned,vector_type_allocator2> * fixation_times,
			  mutation_lookup_table * lookup,
			  const unsigned & generation,const unsigned & twoN )
  {
    static_assert( std::is_base_of<mutation_base,mutation_type>::value,
                   "mutation_type must be derived from KTfwd::mutation_base" );
    typename list_type<mutation_type,list_type_allocator>::iterator i = mutations->begin(),temp;
    while(i != mutations->end())
      {
	assert(i->n <= twoN);			
	i->checked = false;
	if(i->n==twoN )
	  {
	    fixations->push_back(*i);
	    fixation_times->push_back(generation);
	  }
	if( i->n == 0 || i->n == twoN )
	  {
	    lookup->erase(lookup->find(i->pos));
	    temp=i;
	    ++i;
	    mutations->erase(temp);
	  }
	else
	  {
	    ++i;
	  }
      }
  }
  
  template<typename iterator_type>
  void adjust_mutation_counts( iterator_type g , const unsigned & n)
  /*! \brief used internally
    \note Will need a specialization if types have other data that need updating
  */
  {
    auto adjuster = [&n](typename iterator_type::value_type::mutation_list_type_iterator & __m) {
      if(!__m->checked)
	{
	  __m->n=n;
	  __m->checked=true;
	}
      else
	{
	  __m->n += n;
	}
    };
    std::for_each(g->mutations.begin(),g->mutations.end(),
		  std::cref(adjuster));
    std::for_each(g->smutations.begin(),g->smutations.end(),
		  std::cref(adjuster));
   }

  /*! \brief Pick a gamete proportional to its frequency
    Pick a gamete proportional to its frequency
   */
  template< typename gamete_type,
	    typename vector_type_allocator,
	    template <typename,typename> class vector_type > 
  typename vector_type<gamete_type,vector_type_allocator>::iterator
  pgam( gsl_rng * r,
	vector_type<gamete_type,vector_type_allocator > * gametes )
  {
#if  defined(HAVE_BOOST_VECTOR) && !defined(USE_STANDARD_CONTAINERS)
    boost::container::vector<double>freqs(gametes->size(),0);
#else
    std::vector<double>freqs(gametes->size(),0);
#endif
    size_t i = 0;
    std::for_each(gametes->cbegin(),gametes->cend(),
		  [&i,&freqs](const gamete_type & __g) {
		    freqs[i++]=__g.n;
		  });
    fwdpp_internal::gsl_ran_discrete_t_ptr lookup(gsl_ran_discrete_preproc(gametes->size(),&freqs[0]));
    auto rv = gametes->begin()+typename decltype(gametes->begin())::difference_type(gsl_ran_discrete(r,lookup.get()));
    return rv;
  }
										  
  template< typename gamete_type,
	    typename vector_type_allocator,
	    typename mutation_removal_policy,
	    template <typename,typename> class vector_type >
  void update_gamete_list( vector_type<gamete_type,vector_type_allocator > * gametes, 
			   const unsigned & twoN,
			   const mutation_removal_policy & mrp)
  {
    gametes->erase(std::remove_if(gametes->begin(),
     				  gametes->end(),
     				  std::bind(n_is_zero(),std::placeholders::_1)),
     		   gametes->end()); 
#ifndef NDEBUG
#endif
    std::for_each( gametes->begin(),gametes->end(),[&mrp](gamete_type & __g) {
	__g.mutations.erase(std::remove_if( __g.mutations.begin(),
					    __g.mutations.end(),std::cref(mrp) ),
			    __g.mutations.end());
	__g.smutations.erase(std::remove_if( __g.smutations.begin(),
					     __g.smutations.end(),std::cref(mrp) ),
			     __g.smutations.end());
      });
  }

  /*! \brief Multinomial sampling of gametes
    Called by KTfwd::sample_diploid
    \param r GSL random number generator
    \param gbegin Iterator to the start of a vector of gametes
    \param ngametes The size of the vector whose beginning is gbegin
    \param efreqs Vector of frequencies from which to sample.  Must sum to 1 (not checked!).  Typically, this would be the expected frequency of gbegin to (gbegin+i) in the next generation
    \param twoN Twice the number of diploids
   */
  template<typename iterator_type>
  void multinomial_sample_gametes(gsl_rng * r,iterator_type gbegin, 
  				  const size_t & ngametes, 
  				  const double * efreqs, 
  				  const unsigned & twoN)
  {
    unsigned n = twoN;
    double sum=1.;
    for(unsigned i=0;i<ngametes;++i)
      {
    	if( (*(efreqs+i)) > 0. )
    	  {
    	    (gbegin+i)->n = gsl_ran_binomial(r,*(efreqs+i)/sum,n);
    	    sum -= *(efreqs+i);
    	  }
    	else
    	  {
    	    (gbegin+i)->n = 0;
    	  }
    	n -= (gbegin+i)->n;
	adjust_mutation_counts( (gbegin+i), (gbegin+i)->n );
      }
  }
}
#endif /* _UTIL_HPP_ */
