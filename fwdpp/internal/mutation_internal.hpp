#ifndef __FWDPP_INTERNAL_MUTATION_HPP__
#define __FWDPP_INTERNAL_MUTATION_HPP__

#include <algorithm>
#include <cassert>
namespace KTfwd {
  namespace fwdpp_internal
  {
    /*!
      Mechanics of adding a new mutation.

      Ensures that it is always added into a vector
      in less-than-sorted order according to position
    */
    template< typename mcont_t,
	      typename gamete_type>
    void add_new_mutation( const std::size_t idx,
			   const mcont_t & mutations,
			   gamete_type & new_gamete )
    {
      assert(idx<mutations.size());
      if(mutations[idx].neutral)
	{
	  assert(std::find(new_gamete.mutations.begin(),
			   new_gamete.mutations.end(),idx)==new_gamete.mutations.end());
	  new_gamete.mutations.emplace(std::upper_bound(new_gamete.mutations.begin(),
							new_gamete.mutations.end(),mutations[idx].pos,
							[&mutations](const double & __value,const std::size_t __mut) noexcept {
							  assert(__mut<mutations.size());
							  return __value < mutations[__mut].pos;}),
				       idx );
	  assert(gamete_is_sorted_n(new_gamete,mutations));
	}
      else
	{
	  assert(std::find(new_gamete.smutations.begin(),
			   new_gamete.smutations.end(),idx)==new_gamete.smutations.end());
	  new_gamete.smutations.emplace(std::upper_bound(new_gamete.smutations.begin(),
							 new_gamete.smutations.end(),mutations[idx].pos,
							 [&mutations](const double & __value,const std::size_t __mut) noexcept {
							   assert(__mut<mutations.size());
							   return __value < mutations[__mut].pos;}),
					idx );
	  assert(gamete_is_sorted_s(new_gamete,mutations));
	}

    }
    
    template<typename mmodel,
    	     typename gamete_type,
    	     typename mlist_type,
    	     typename queue_t>
    inline typename std::result_of<mmodel(queue_t &,mlist_type &)>::type
    mmodel_dispatcher( const mmodel & m, gamete_type & , mlist_type & mutations, queue_t & recycling_bin)
    {
      return m(recycling_bin,mutations);
    }

    template<typename mmodel,
	     typename gamete_type,
    	     typename mlist_type,
    	     typename queue_t>
    inline typename std::result_of<mmodel(queue_t &,gamete_type &,mlist_type &)>::type
    mmodel_dispatcher( const mmodel & m, gamete_type & g, mlist_type & mutations, queue_t & recycling_bin)
    {
      return m(recycling_bin,g,mutations);
    }

    /*!
      Apply mutation model N times to a new gamete.
      Updates mutation list
    */
    template<typename queue_type,
	     typename mutation_model,
    	     typename mlist_type,
    	     typename gamete_type>
    void add_N_mutations_recycle( queue_type & recycling_bin,
				  const mutation_model & mmodel,
				  const unsigned & n,
				  mlist_type & mutations,
				  gamete_type & g)
    {
      assert(gamete_is_sorted_n(g,mutations));
      assert(gamete_is_sorted_s(g,mutations));
      for( unsigned i = 0 ; i < n ; ++i )
    	{
	  add_new_mutation(mmodel_dispatcher(mmodel,g,mutations,recycling_bin),mutations,g);
    	}
    }
  }
}

#endif
