#ifndef FWDPP_INTERNAL_RECYCLING
#define FWDPP_INTERNAL_RECYCLING

#include <queue>

namespace KTfwd
{
  namespace fwdpp_internal
  {
    template<class T> using recycling_bin_t = std::queue<T>;
    
    template<typename mlist_t>
    recycling_bin_t<typename mlist_t::iterator> make_mut_queue( mlist_t * mutations )
    {
      recycling_bin_t<typename mlist_t::iterator> rv;
      for(auto mitr = mutations->begin();mitr!=mutations->end();++mitr)
	{
	  if(!mitr->n && !mitr->checked) rv.push(mitr);
	}
      return rv;
    }

    template<typename glist_t>
    recycling_bin_t<typename glist_t::iterator> make_gamete_queue( glist_t * gametes )
    {
      recycling_bin_t<typename glist_t::iterator> rv;
      for(auto gitr = gametes->begin();gitr!=gametes->end();++gitr)
	{
	  if(!gitr->n) rv.push(gitr);
	}
      return rv;
    }

    template<typename iterator_type,
	     typename glist_t,
	     typename queue_t,
	     typename lookup_type>
    inline void recycle_gamete( iterator_type & g1,
				glist_t * gametes,
				queue_t & gamete_recycling_bin, lookup_type & gamete_lookup,
				typename iterator_type::value_type::mutation_container & neutral,
				typename iterator_type::value_type::mutation_container & selected )
    {
      //Try to recycle
      if( ! gamete_recycling_bin.empty() )
	{
	  g1 = gamete_recycling_bin.front();
	  assert(!g1->n);
	  g1->n=0u;
	  std::swap(g1->mutations,neutral);
	  std::swap(g1->smutations,selected);
	  gamete_recycling_bin.pop();
	}
      else
	{
	  g1 = gametes->emplace(gametes->end(),0u,std::move(neutral),std::move(selected));
	}
      gamete_lookup.update(g1);
    }

    /*!
      \brief Helper function for mutation policies

      This function minimizes code duplication when writing mutation models.  It abstracts
      the operations needed to recycle an extinct mutation.

      \param mutation_recycling_bin  A FIFO queue of iterators pointing to extinct mutations.
      \param mutations A list of mutation objects
      \param args Parameter pack to be passed to constructor of an mlist_t::value_type
     */
    template<typename queue_t,
	     typename mlist_t,
	     class... Args >
    typename mlist_t::iterator recycle_mutation_helper( queue_t & mutation_recycling_bin,
							mlist_t * mutations,
							Args&&... args )
    {
      if(!mutation_recycling_bin.empty())
	{
	  auto rv = mutation_recycling_bin.front();
	  mutation_recycling_bin.pop();
	  *rv = typename mlist_t::value_type(args...);
	  return rv;
	}
      return mutations->emplace(mutations->end(),std::forward<Args>(args)...);
    }
  }
}

#endif
