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
      
      fwdpp 0.3.6 further changed the lookup scheme to have a
      friendlier API as well as a much faster lookup for the case of 
      simulations involving selected mutations
    */
    template<typename gcont_t>
    struct gamete_lookup {
      using gcont_t_itr = typename gcont_t::iterator;
      using lookup_table_t = std::multimap<double,gcont_t_itr>;
      using result_type = std::pair<typename lookup_table_t::iterator,typename lookup_table_t::iterator>;
      using inner_t = typename lookup_table_t::value_type;
      lookup_table_t lookup_table;

      inline void update_details( gcont_t_itr g ) 
      {
	double npos0 = ((g->mutations.empty()) ? -std::numeric_limits<double>::max() : g->mutations[0]->pos);
	double spos0 = ((g->smutations.empty()) ? -std::numeric_limits<double>::max() : g->smutations[0]->pos);
	lookup_table.insert( std::make_pair( npos0*double(g->mutations.size()) + spos0*double(g->smutations.size()), g ) );
      }

      explicit gamete_lookup(gcont_t * gametes) : lookup_table( lookup_table_t() )
      {
	for( auto g = gametes->begin() ; g != gametes->end() ; ++g )
	  {
	    update_details(g);
	  }
      }

      inline result_type lookup( const typename gcont_t::value_type & g ) 
      {
	double npos0 = ((g.mutations.empty()) ? -std::numeric_limits<double>::max() : g.mutations[0]->pos);
	double spos0 = ((g.smutations.empty()) ? -std::numeric_limits<double>::max() : g.smutations[0]->pos);
	return lookup_table.equal_range(  npos0*double(g.mutations.size()) + spos0*double(g.smutations.size()) );
      }

      inline void update( gcont_t_itr g ) 
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
