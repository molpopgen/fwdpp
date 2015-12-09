#ifndef FWDPP_INTERNAL_SAMPLE_DIPLOID_HELPERS
#define FWDPP_INTERNAL_SAMPLE_DIPLOID_HELPERS

namespace KTfwd
{
  namespace fwdpp_internal
  {
    template<typename glist_t>
    inline void process_glist( glist_t * gametes )
    {
      for(auto itr = gametes->begin() ; itr != gametes->end() ; ++itr )
	{
	  if(itr->n) adjust_mutation_counts(itr,itr->n);
	}
    }
  }
}

#endif
