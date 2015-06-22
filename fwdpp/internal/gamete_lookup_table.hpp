#ifndef __FWDPP__INTERNAL_GAMETE_LOOKUP_TABLE_HPP__
#define __FWDPP__INTERNAL_GAMETE_LOOKUP_TABLE_HPP__

#include <vector>
#include <map>
#include <unordered_map>

namespace KTfwd {
  namespace fwdpp_internal {
#ifdef VECTOR_GLOOKUP
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
#elif defined MULTIMAP_GLOOKUP
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
#else //unordered multimap
    template< typename gcont_t >
    std::unordered_multimap<std::uint32_t,typename gcont_t::iterator>
    gamete_lookup_table(gcont_t * gametes )
    {
      using rv_t = std::unordered_multimap<std::uint32_t,typename gcont_t::iterator>;
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
