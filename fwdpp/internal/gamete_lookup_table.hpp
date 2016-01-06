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
    template<typename gcont_t,
	     typename mcont_t>
    struct gamete_lookup
    {
      using lookup_table_t = std::multimap<double,size_t>;
      using result_type = std::pair< lookup_table_t::iterator, lookup_table_t::iterator>;
      using inner_t = typename lookup_table_t::value_type;
      lookup_table_t lookup_table;

      inline double keyit( const std::vector<size_t> & mc,
			   const mcont_t & mutations ) const
      {
	return (mc.empty()) ? -std::numeric_limits<double>::max() : mutations[mc[0]].pos;
      }

      inline void update_details( size_t g,
				  const gcont_t & gametes,
				  const mcont_t & mutations)
      {
	lookup_table.emplace( std::make_pair( keyit(gametes[g].mutations,mutations)*double(gametes[g].mutations.size()) +
					      keyit(gametes[g].smutations,mutations)*double(gametes[g].smutations.size()), g) );
      }

      inline result_type lookup( const std::vector<size_t> & n,
				 const std::vector<size_t> & s,
				 const mcont_t & mutations )
      {
	return lookup_table.equal_range(  keyit(n,mutations)*double(n.size()) + keyit(s,mutations)*double(s.size()) );
      }

      inline void update(size_t idx,
			 const gcont_t & gametes,
			 const mcont_t & mutations)
      {
	update_details(idx,gametes,mutations);
      };

      explicit gamete_lookup( const gcont_t & gametes,
			      const mcont_t & mutations )
      {
	for(size_t g=0;g<gametes.size();++g)
	  {
	    if(gametes[g].n)
	      {
		update_details(g,gametes,mutations);
	      }
	  }
      }
    };

    /*!
      Convenience function to return a lookup table.
      Exists so that calling environ can use auto and not
      deal w/template type, which is an attempt at future-proofing.
    */
    template< typename gcont_t, typename mcont_t >
    gamete_lookup<gcont_t,mcont_t>
    gamete_lookup_table(const gcont_t & gametes,
			const mcont_t & mutations)  {
      return gamete_lookup<gcont_t,mcont_t>(gametes,mutations);
    }
  }
}

#endif
