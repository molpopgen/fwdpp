/// \file util.hpp
#ifndef _UTIL_HPP_
#define _UTIL_HPP_

#include <fwdpp/forward_types.hpp>
#include <fwdpp/fwd_functional.hpp>
#include <fwdpp/type_traits.hpp>
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
    static_assert( typename traits::is_mutation_t<mutation_type>::type(),
                   "mutation_type must be derived from KTfwd::mutation_base" );
    std::for_each(mutations->begin(),mutations->end(),[](mutation_type & __m){__m.checked=false;});
  }

  /*
    Label all extinct variants for recycling

    \note: lookup must be compatible with lookup->erase(lookup->find(double))
  */
  template<typename mutation_list_type,
	   typename mutation_lookup_table>
  void update_mutations( mutation_list_type * mutations, 
			 mutation_lookup_table * lookup )
  {
    static_assert( typename traits::is_mutation_t<typename mutation_list_type::value_type>::type(),
		   "mutation_type must be derived from KTfwd::mutation_base" );
    for(auto i=mutations->begin();i!=mutations->end();++i)
      {
	if(!i->checked)
	  {
	    i->n=0;
	    lookup->erase(i->pos);
	  }
	i->checked=false;
      }
  }

  /*!
    Label all extinct and fixed variants for recycling

    \note: lookup must be compatible with lookup->erase(lookup->find(double))
  */
  template<typename mutation_list_type,
	   typename mutation_lookup_table>
  void update_mutations( mutation_list_type & mutations, 
			 mutation_lookup_table & lookup,
			 const std::vector<uint_t> & mcounts,
			 const unsigned twoN)
  {
    static_assert( typename traits::is_mutation_t<typename mutation_list_type::value_type>::type(),
		   "mutation_type must be derived from KTfwd::mutation_base" );
    for(std::size_t i = 0 ; i < mcounts.size() ; ++i)
      {
	assert(mcounts[i] <= twoN);
	if(mcounts[i]==twoN || !mcounts[i] )
	  {
	    lookup.erase(mutations[i].pos);
	    mcounts[i]=0;
	  }
      }
  }

  /*!
    Label all extinct variants for recycling, copy fixations and fixation times
    into containers.

    \note: lookup must be compatible with lookup->erase(lookup->find(double))
  */
  template<//typename mutation_type,
	   typename mutation_list_type,
	   typename fixation_container_t,
	   typename fixation_time_container_t,
	   //typename vector_type_allocator1,
	   //typename vector_type_allocator2,
	   //typename list_type_allocator,
	   //template <typename,typename> class vector_type,
	   //template <typename,typename> class list_type,
	   typename mutation_lookup_table>
  void update_mutations( mutation_list_type & mutations, 
			 fixation_container_t & fixations, 
			 fixation_time_container_t & fixation_times,
			 mutation_lookup_table & lookup,
			 const std::vector<uint_t> & mcounts,
			 const unsigned & generation,const unsigned & twoN )
  {
    static_assert( typename traits::is_mutation_t<typename mutation_list_type::value_type>::type(),
		   "mutation_type must be derived from KTfwd::mutation_base" );
    for(unsigned i=0;i<mcounts.size();++i)
      {
	assert(mcounts[i] <= twoN);
	if(mcounts[i]==twoN )
	  {
	    fixations.push_back(mutations[i]);
	    fixation_times.push_back(generation);
	    lookup.erase(mutations[i].pos);
	  }
      }
  }
  
  /*!
    Add elements to a container via emplace_back
    
    \param c The container
    \param N final desired size of the container
    \param args Arguments for the constructor of container_type::value_type

    \note if N <= c.size(), nothing happens.
  */
  template<typename container_type,
	   class... Args>
  void add_elements( container_type & c,
		     size_t N,
		     Args&&... args )
  {
    for( size_t i = c.size() ; i < N ; ++i )
      {
	c.emplace_back(std::forward<Args>(args)...);
      }
  }
}
#endif /* _UTIL_HPP_ */
