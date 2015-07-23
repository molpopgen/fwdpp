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

namespace KTfwd
{
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
  */
  template<typename mutation_type,
	   typename list_type_allocator,
	   template <typename,typename> class list_type >
  void remove_lost( list_type<mutation_type,list_type_allocator> * mutations )
  {
    static_assert( std::is_base_of<mutation_base,mutation_type>::value,
                   "mutation_type must be derived from KTfwd::mutation_base" );
    for(auto i = mutations->begin(); i != mutations->end() ; )
      {
	if(!i->checked)
	  {
	    i = mutations->erase(i);
	  }
	  else
	  {
	    i->checked = false;
	    ++i;
	  }
      }
  }

 /*! \brief Remove mutations from population
    Removes mutations that are lost.
    \note: lookup must be compatible with lookup->erase(lookup->find(double))
   */
  template<typename mutation_type,
	   typename list_type_allocator,
	   template <typename,typename> class list_type,
	   typename mutation_lookup_table>
  void remove_lost( list_type<mutation_type,list_type_allocator> * mutations, mutation_lookup_table * lookup )
  {
    static_assert( std::is_base_of<mutation_base,mutation_type>::value,
                   "mutation_type must be derived from KTfwd::mutation_base" );
    for(auto i = mutations->begin() ; i != mutations->end() ; )
      {
	if(!i->checked)
	  {
	    lookup->erase(lookup->find(i->pos));
	    i=mutations->erase(i);
	  }
	else
	  {
	    i->checked=false;
	    ++i;
	  }
      }
  }

  /*! \brief Remove mutations from population
    Removes mutations that are fixed or lost.
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
    for(auto i = mutations->begin() ; i != mutations->end() ; )
      {
	assert(i->n <= twoN);			
	if(i->n==twoN )
	  {
	    fixations->push_back(*i);
	    fixation_times->push_back(generation);
	  }
	if( !i->checked || i->n == twoN )
	  {
	    i=mutations->erase(i);
	  }
	else
	  {
	    i->checked=false;
	    ++i;
	  }
      }
  }

  /*! \brief Remove mutations from population
    Removes mutations that are fixed or lost.
    \note: lookup must be compatible with lookup->erase(lookup->find(double))
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
    for(auto i=mutations->begin();i!=mutations->end();)
      {
	assert(i->n <= twoN);
	if(i->n==twoN )
	  {
	    fixations->push_back(*i);
	    fixation_times->push_back(generation);
	  }
	if(!i->checked ||  i->n == twoN )
	  {
	    lookup->erase(lookup->find(i->pos));
	    i = mutations->erase(i);
	  }
	else
	  {
	    i->checked=false;
	    ++i;
	  }
      }
    
  }

  template<typename iterator_type>
  void adjust_mutation_counts( iterator_type & g , const unsigned & n)
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
}
#endif /* _UTIL_HPP_ */
