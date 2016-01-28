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
      //Establish a pointer to the container to be modified
      auto mc = (mutations[idx].neutral) ? &new_gamete.mutations : &new_gamete.smutations;

      assert( std::find(mc->cbegin(),mc->cend(),idx) == mc->end() );
      //Insert new mutation at position retaining sort-by-position order
      mc->emplace(std::upper_bound(mc->begin(),
				   mc->end(),mutations[idx].pos,
				   [&mutations,&idx](const double & __value,const std::size_t __mut) noexcept {
				     assert(__mut<mutations.size());
				     return __value < mutations[__mut].pos;}),
		  idx );

      assert(std::is_sorted(mc->cbegin(),mc->cend(),
			    [&mutations](const std::size_t i, const std::size_t j)noexcept
			    {
			      return mutations[i].pos<mutations[j].pos;
			    }));
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
