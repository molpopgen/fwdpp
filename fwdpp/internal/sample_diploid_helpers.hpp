#ifndef FWDPP_INTERNAL_SAMPLE_DIPLOID_HELPERS
#define FWDPP_INTERNAL_SAMPLE_DIPLOID_HELPERS

#include <algorithm>

namespace KTfwd
{
  namespace fwdpp_internal
  {
    template<typename iterator_type>
    void adjust_mutation_counts( iterator_type & g , const unsigned & n)
    /*! 
      Update the counts of each mutation.

      This function is called by KTfwd::fwdpp_internal::adjust_mutation_counts
    */
    {
      auto adjuster = [&n](typename iterator_type::value_type::mutation_list_type_iterator & __m) {
	if(!__m->checked)
	  {
	    __m->n=n;
	    __m->checked=true;
	  }
	else
	  {
	    __m->n += n;
	  }
      };
      std::for_each(g->mutations.begin(),g->mutations.end(),
		    std::cref(adjuster));
      std::for_each(g->smutations.begin(),g->smutations.end(),
		    std::cref(adjuster));
    }
    
    template<typename glist_t>
    inline void process_glist( glist_t * gametes )
    /*!
      For every non-extinct gamete, increment the counts of its mutations
      via a call to KTfwd::fwdpp_internal::adjust_mutation_counts.
    */
    {
      for(auto itr = gametes->begin() ; itr != gametes->end() ; ++itr )
	{
	  if(itr->n) adjust_mutation_counts(itr,itr->n);
	}
    }
  }
}

#endif
