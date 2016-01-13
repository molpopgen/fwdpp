#ifndef __FWDPP_INTERNAL_GAMETE_CLEANER_HPP__
#define __FWDPP_INTERNAL_GAMETE_CLEANER_HPP__

#include <type_traits>

namespace KTfwd {
  namespace fwdpp_internal {
    /*! \brief Handles removal of pointers to mutations from gametes after sampling
      Intended use is when std::is_same< mutation_removal_policy, KTfwd::remove_nothing >::type is true.
      Called by KTfwd::sample_diploid via dispatch.
    */
    template<typename gcont_t, typename mutation_removal_policy>
    inline void gamete_cleaner(gcont_t &,const std::vector<uint_t> &, const mutation_removal_policy &, std::true_type)
    {
      return;
    }

    /*! \brief Handles removal of pointers to mutations from gametes after sampling
      Intended use is when std::is_same< mutation_removal_policy, KTfwd::remove_nothing >::type is false.
      Called by KTfwd::sample_diploid via dispatch.
    */
    template<typename gcont_t>//, typename mutation_removal_policy>
    inline void gamete_cleaner(gcont_t & gametes, const std::vector<uint_t> & mcounts,
			       const uint_t twoN, std::false_type)
    {
      for( auto & g : gametes )
	{
	  if(g.n)
	    {
	      g.mutations.erase(std::remove_if(g.mutations.begin(),
					       g.mutations.end(),
					       [&mcounts,&twoN](const std::size_t & i) noexcept
					       {
						 return mcounts[i]==twoN;
					       }),
				g.mutations.end());
	      g.smutations.erase(std::remove_if(g.smutations.begin(),
						g.smutations.end(),
						[&mcounts,&twoN](const std::size_t & i) noexcept
						{
						  return mcounts[i]==twoN;
						}),
				 g.smutations.end());
	    }
	}
    }
  }
}

#endif
