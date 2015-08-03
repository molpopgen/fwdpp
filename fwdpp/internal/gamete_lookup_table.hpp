#ifndef __FWDPP__INTERNAL_GAMETE_LOOKUP_TABLE_HPP__
#define __FWDPP__INTERNAL_GAMETE_LOOKUP_TABLE_HPP__

#include <map>
#include <limits>
namespace KTfwd {
  namespace fwdpp_internal {

    /*!
      Allows KTfwd::recombine_gametes to quickly
      decide if a recombinant is a gamete that does, 
      or does not, currently exist in the population.

      Introduced in fwdpp 0.3.5, replacing the multimap
      introduced in 0.3.3 that simply associated iterators
      to gametes with the number of mutations.
      
      The new lookup separates out the number of neutral
      and selected mutations, which speeds up simulations
      with selection.
    */
    template<typename gcont_t>
    struct gamete_lookup {
      using uint_t = std::uint32_t;
      using gcont_t_itr = typename gcont_t::iterator;
      using mmap_t = std::multimap<std::pair<double,double>,gcont_t_itr>;
      using inner_t = typename mmap_t::value_type;
      using lookup_table_t = std::map<unsigned,
				      std::map<unsigned,
					       mmap_t
					       >
				      >;

      using map_t = std::map<unsigned,mmap_t>;
      using result_type = std::pair<bool,std::pair<typename mmap_t::iterator,typename mmap_t::iterator> >;
      lookup_table_t lookup_table;
      
      void update_details( gcont_t_itr g ) 
      {
	double npos0 = ((g->mutations.empty()) ? -std::numeric_limits<double>::max() : g->mutations[0]->pos);
	double spos0 = ((g->smutations.empty()) ? -std::numeric_limits<double>::max() : g->smutations[0]->pos);
	const auto pospair = std::make_pair(spos0,npos0);
	auto __outer = lookup_table.find(g->smutations.size());
	 if( __outer == lookup_table.end() ) 
	   {
	     lookup_table[g->smutations.size()] = map_t({std::make_pair(g->mutations.size(),mmap_t({std::make_pair(pospair,g)}))});
	   }
	 else {
	   auto __i = __outer->second.find( g->mutations.size() );
	   if(__i == __outer->second.end())
	     {
	       __outer->second.insert( std::make_pair(g->mutations.size(),mmap_t({std::make_pair(pospair,g)})) );
	     }
	   else {
	     __i->second.insert( std::make_pair(pospair,g) );
	   }
	 }
      }

      explicit gamete_lookup(gcont_t * gametes) : lookup_table( lookup_table_t() )
      {
	for( auto g = gametes->begin() ; g != gametes->end() ; ++g )
	  {
	    update_details(g);
	  }
      }

      result_type lookup( const typename gcont_t::value_type & g ) 
      {
	auto itr = lookup_table.find(g.smutations.size());
	if(itr == lookup_table.end())
	  {
	    return std::make_pair(false,std::make_pair(typename mmap_t::iterator(),typename mmap_t::iterator()));
	  }
	auto jtr = itr->second.find(g.mutations.size());
	if( jtr == itr->second.end() )
	  {
	    return std::make_pair(false,std::make_pair(typename mmap_t::iterator(),typename mmap_t::iterator()));
	  }
	double npos0 = ((g.mutations.empty()) ? -std::numeric_limits<double>::max() : g.mutations[0]->pos);
	double spos0 = ((g.smutations.empty()) ? -std::numeric_limits<double>::max() : g.smutations[0]->pos);
	return std::make_pair(true,jtr->second.equal_range(std::make_pair(spos0,npos0)));
      }

      void update( gcont_t_itr g ) 
      {
	update_details(g);
      }
    };

    /*!
      Convenience function to return a lookup table.
      Exists so that calling environ can use auto and not 
      deal w/template type, which is an attempt at future-proofing.
    */
    template< typename gcont_t >
    gamete_lookup<gcont_t>
    gamete_lookup_table(gcont_t * gametes )  {
      return gamete_lookup<gcont_t>(gametes);
    }
  }
}

#endif
