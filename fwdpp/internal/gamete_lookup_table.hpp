#ifndef __FWDPP__INTERNAL_GAMETE_LOOKUP_TABLE_HPP__
#define __FWDPP__INTERNAL_GAMETE_LOOKUP_TABLE_HPP__

#ifdef VECTOR_GLOOKUP
#include <vector>
#else
#include <map>
#endif

namespace KTfwd {
  namespace fwdpp_internal {
#ifdef FWDPP_VECTOR_GLOOKUP
    template< typename gcont_t >
    std::vector<std::pair<std::uint32_t,typename gcont_t::iterator> >
    gamete_lookup_table(gcont_t * gametes )
    {
      using pair_t = std::pair<std::uint32_t,typename gcont_t::iterator>;
      using rv_t = std::vector<pair_t>;
      rv_t rv;
      rv.reserve(2*gametes->size());
      for(auto itr = gametes->begin() ; itr != gametes->end() ; ++itr )
	{
	  rv.emplace_back(std::make_pair(itr->mutations.size()+itr->smutations.size(),itr) );
	}
      std::sort(rv.begin(),rv.end(),[](const pair_t & a,const pair_t & b) { return a.first<b.first; });
      return rv;
    }
#else
    template< typename gcont_t >
    std::multimap<std::uint32_t,typename gcont_t::iterator>
    gamete_lookup_table(gcont_t * gametes )
    {
      using rv_t = std::multimap<std::uint32_t,typename gcont_t::iterator>;
      rv_t rv;
      for(auto itr = gametes->begin() ; itr != gametes->end() ; ++itr )
	{
	  rv.insert( std::make_pair(itr->mutations.size()+itr->smutations.size(),itr) );
	}
      return rv;
    }
#endif
  }
}

#endif
