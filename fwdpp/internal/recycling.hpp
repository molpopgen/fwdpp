#ifndef FWDPP_INTERNAL_RECYCLING
#define FWDPP_INTERNAL_RECYCLING

#include <queue>

namespace KTfwd
{
  namespace fwdpp_internal
  {
    template<typename mlist_t>
    typename std::queue<typename mlist_t::iterator> make_mut_queue( mlist_t * mutations )
    {
      std::queue<typename mlist_t::iterator> rv;
      for(auto mitr = mutations->begin();mitr!=mutations->end();++mitr)
	{
	  if(!mitr->n && !mitr->checked) rv.push(mitr);
	}
      return rv;
    }

    template<typename glist_t>
    typename std::queue<typename glist_t::iterator> make_gamete_queue( glist_t * gametes )
    {
      std::queue<typename glist_t::iterator> rv;
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

    template<typename queue_t,
	     typename mlist_t,
	     class ... Args >
    typename mlist_t::iterator mutation_helper( queue_t & mutation_recycling_bin,
						mlist_t * mutations,
						Args ... args )
    {
      if(!mutation_recycling_bin.empty())
	{
	  auto rv = mutation_recycling_bin.front();
	  mutation_recycling_bin.pop();
	  *rv = typename mlist_t::value_type(args...);
	  return rv;
	}
      return mutations->emplace(mutations->end(),args...);
    }
  }
}

#endif
