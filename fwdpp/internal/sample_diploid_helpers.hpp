#ifndef FWDPP_INTERNAL_SAMPLE_DIPLOID_HELPERS
#define FWDPP_INTERNAL_SAMPLE_DIPLOID_HELPERS

#include <algorithm>

namespace KTfwd
{
  namespace fwdpp_internal
  {
    template<typename gcont_t,
	     typename mcont_t>
    inline void process_gametes( const gcont_t & gametes,
				 const mcont_t & mutations,
				 std::vector<uint_t> & mcounts)
    /*!
      For every non-extinct gamete, increment the counts of its mutations
      via a call to KTfwd::fwdpp_internal::adjust_mutation_counts.
    */
    {
      if(mutations.size()>mcounts.size())
	{
	  mcounts.resize(mutations.size(),0);
	}
      //zero out mcounts
      for(auto & mc : mcounts) mc=0;
      //update mutation counts
      for(const auto & g : gametes)
	{
	  if(g.n) //only do this for extant gametes
	    {
	      for(const auto & m : g.mutations) mcounts[m]+=g.n;
	      for(const auto & m : g.smutations) mcounts[m]+=g.n;
	    }
	}
    }
  }
}

#endif
