#ifndef __FWDPP_INTERNAL_GAMETE_CLEANER_HPP__
#define __FWDPP_INTERNAL_GAMETE_CLEANER_HPP__

#include <type_traits>
#include <fwdpp/fwd_functional.hpp>
namespace KTfwd {
  namespace fwdpp_internal {
    /*! \brief Handles removal of pointers to mutations from gametes after sampling
      Intended use is when std::is_same< mutation_removal_policy, KTfwd::remove_nothing >::type is true.
      Called by KTfwd::sample_diploid
    */
    template<typename gcont_t,typename mcont_t,typename mutation_removal_policy>
    inline
    typename std::enable_if< std::is_same<mutation_removal_policy,KTfwd::remove_nothing>::value >::type
    gamete_cleaner(gcont_t &,const mcont_t &,const std::vector<uint_t> &, const uint_t, const mutation_removal_policy &)
    {
      return;
    }

    /*! \brief Handles removal of pointers to mutations from gametes after sampling
      Intended use is when std::is_same< mutation_removal_policy, KTfwd::true_type >::type is true.
      Called by KTfwd::sample_diploid
    */
    template<typename gcont_t,typename mcont_t,typename mutation_removal_policy>
    inline typename std::enable_if< std::is_same<mutation_removal_policy,std::true_type>::value >::type
    gamete_cleaner(gcont_t & gametes, const mcont_t &, const std::vector<uint_t> & mcounts,
		   const uint_t twoN, const mutation_removal_policy &)
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

    /*! \brief Handles removal of pointers to mutations from gametes after sampling
      This overload handles truly custom policies, which must take a mutation type as an argument.
      Called by KTfwd::sample_diploid
    */
    template<typename gcont_t,typename mcont_t,typename mutation_removal_policy>
    inline typename std::enable_if<
      !std::is_same<mutation_removal_policy,std::true_type>::value &&
    !std::is_same<mutation_removal_policy,KTfwd::remove_nothing>::value
    >::type
    gamete_cleaner(gcont_t & gametes, const mcont_t & mutations,const std::vector<uint_t> & mcounts,
		   const uint_t twoN, const mutation_removal_policy & mp)
    {
      for( auto & g : gametes )
    	{
    	  if(g.n)
    	    {
    	      g.mutations.erase(std::remove_if(g.mutations.begin(),
    					       g.mutations.end(),
    					       [&mutations,&mcounts,&twoN,&mp](const std::size_t & i) noexcept
    					       {
    						 return mcounts[i]==twoN && mp(mutations[i]);
    					       }),
    				g.mutations.end());
    	      g.smutations.erase(std::remove_if(g.smutations.begin(),
    						g.smutations.end(),
    						[&mutations,&mcounts,&twoN,&mp](const std::size_t & i) noexcept
    						{
    						  return mcounts[i]==twoN && mp(mutations[i]);
    						}),
    				 g.smutations.end());
    	    }
    	}
    }
  }
}

#endif
