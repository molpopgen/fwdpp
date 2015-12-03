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
  }
}

#endif
