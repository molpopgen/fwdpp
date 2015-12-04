#ifndef __FWDPP_INTERNAL_MUTATION_HPP__
#define __FWDPP_INTERNAL_MUTATION_HPP__

#include <algorithm>
#include <type_traits>

namespace KTfwd {
  namespace fwdpp_internal
  {
    /*!
      Mechanics of adding a new mutation.

      Ensures that it is always added into a vector
      in less-than-sorted order according to position
    */
    template< typename gamete_type,
	      typename mutation_iterator>
    void add_new_mutation( mutation_iterator mitr,
			   gamete_type & new_gamete )
    {
      if(mitr->neutral)
	{
	  new_gamete.mutations.emplace(std::upper_bound(new_gamete.mutations.begin(),
							new_gamete.mutations.end(),mitr->pos,
							[](const double & __value,const mutation_iterator & __mut){
							  return __value < __mut->pos;}),
				       mitr );
	}
      else
	{
	  new_gamete.smutations.emplace(std::upper_bound(new_gamete.smutations.begin(),
							 new_gamete.smutations.end(),mitr->pos,
							 [](const double & __value,const mutation_iterator & __mut){return __value < __mut->pos;}),
					mitr );
	}
    }
    
    template<typename mmodel,
    	     typename gamete_type,
    	     typename mlist_type,
    	     typename queue_t>
    inline
    typename std::enable_if< std::is_same< typename std::result_of<mmodel(queue_t &,mlist_type *)>::type , typename mlist_type::iterator >::value,
     			     typename mlist_type::iterator >::type
    mmodel_dispatcher( const mmodel & m, gamete_type & g, mlist_type * mutations, queue_t & recycling_bin)
    {
      return m(recycling_bin,mutations);
    }

    template<typename mmodel,
	     typename gamete_type,
    	     typename mlist_type,
    	     typename queue_t>
    inline
    typename std::enable_if< std::is_same< typename std::result_of<mmodel(queue_t &,gamete_type &,mlist_type *)>::type , typename mlist_type::iterator >::value,
     			     typename mlist_type::iterator >::type
    mmodel_dispatcher( const mmodel & m, gamete_type & g, mlist_type * mutations, queue_t & recycling_bin)
    {
      return m(recycling_bin,g,mutations);
    }
    
    /*!
      Apply mutation model N times to a new gamete.
      Updates mutation list
    */
    template<typename mutation_model,
    	     typename mutation_insertion_policy,
    	     typename mlist_type,
    	     typename gamete_type>
    void add_N_mutations( const mutation_model & mmodel,
    			  const mutation_insertion_policy & mpolicy,
    			  const unsigned & n,
    			  mlist_type * mutations,
    			  gamete_type & g)
    {
      for( unsigned i = 0 ; i < n ; ++i )
    	{
	  auto m = mmodel_dispatcher(mmodel,g,mutations);
    	  auto mitr = mpolicy(std::move(m),mutations);
    	  add_new_mutation(mitr,g);
    	}
    }

    /*!
      Apply mutation model N times to a new gamete.
      Updates mutation list
    */
    template<typename queue_type,
	     typename mutation_model,
    	     typename mutation_insertion_policy,
    	     typename mlist_type,
    	     typename gamete_type>
    void add_N_mutations_recycle( queue_type & recycling_bin,
				  const mutation_model & mmodel,
				  const mutation_insertion_policy & mpolicy,
				  const unsigned & n,
				  mlist_type * mutations,
				  gamete_type & g)
    {
      for( unsigned i = 0 ; i < n ; ++i )
    	{
	  add_new_mutation(mmodel_dispatcher(mmodel,g,mutations,recycling_bin),g);
    	}
    }
  }
}

#endif
