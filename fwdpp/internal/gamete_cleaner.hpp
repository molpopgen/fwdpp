#ifndef __FWDPP_INTERNAL_GAMETE_CLEANER_HPP__
#define __FWDPP_INTERNAL_GAMETE_CLEANER_HPP__

#include <vector>
#include <algorithm>
#include <type_traits>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/fwd_functional.hpp>

/*!
  \file gamete_cleaner.hpp

  This file handles the "pruning" of fixations from gametes.

  Whether or not any pruning happens depends on policies.

  Version 0.4.7 of the library made a major change resulting in improved 
  run-time performance, especially for low mutation rates.

  Consider the case of neutral mutations occurring at rate mu (per site,
  per generation).  The rate of fixation is also mu.  Thus, when mu is low,
  the rate of fixation is low, and the expected time in generations between 
  fixations is 1/mu.  In other words, calling gamete_cleaner every generation
  is an unnecessary run-time cost.  Now, we check if any fixations exist, and return
  if there are none.  We check separately for neutral- and non-neutral fixations.

  Naively, we could simply scan the mutations container of the entire simulation for fixations.
  But, we can do better than that because:

  1. Fixations are, by definition, found in all gametes.  Thus, it suffices to check the
  first (extant) gamete.

  2. When we find the first fixation in a gamete, it is by definition the fixation with the smallest
  position, out of all fixations.  We can store its value and use it to look up the first fixation
  in all of the remaining gametes.  Searching for the first fixation in this way avoids cache misses that 
  are unavoidable when we do out-of-order lookups in "mcounts" for the remaining fixations.
*/

namespace KTfwd {
  namespace fwdpp_internal {

    //First, we have a set of helper functions:

    //! Wrapper around std::find_if to find next gamete where n > 0
    template<typename gcont_t_iterator>
    inline gcont_t_iterator next_extant_gamete(gcont_t_iterator beg,
					       gcont_t_iterator end) noexcept
    {
      return std::find_if(beg,end,[](const typename gcont_t_iterator::value_type & g) noexcept { return g.n; });
    }

    //! Find the first mutation in a gamete that is also a fixation
    template<typename mut_index_cont,typename mcounts_t>
    inline typename mut_index_cont::const_iterator find_fixation(const mut_index_cont & mc,
								 const mcounts_t & mcounts,
								 const uint_t twoN) noexcept
    {
      return std::find_if(mc.cbegin(),mc.cend(),[&mcounts,&twoN](const std::size_t & i) noexcept
			  {
			    return mcounts[i]==twoN;
			  });
    }

    //! Wrapper for erase/remove idiom.
    template<typename mut_index_cont,typename mcounts_t>
    inline void gamete_cleaner_erase_remove_idiom_wrapper( mut_index_cont & mc,
							   const mcounts_t & mcounts,
							   const typename mut_index_cont::value_type & first_fixation,
							   const uint_t twoN ) noexcept
    {
      /*
	The first call to std::find relies on a de-referencing of the return value of 
	find_fixation, passed to this function as first_fixation.  Because mc is sorted according to mutation position, 
	first_fixation is the fixation with the smallest position.  Using std::find like this avoids some calls to the lambda
	expression, which experiences cache-misses because mcounts is NOT sorted according to mutation position.
      */
      mc.erase(std::remove_if(std::find(mc.begin(),mc.end(),first_fixation),
			      mc.end(),[&mcounts,&twoN](const typename mut_index_cont::value_type & i) noexcept
			      {
				return mcounts[i]==twoN;
			      }),mc.end());
    }

    //! Find the first mutation in a gamete that is also a fixation AND satisfies the policy
    template<typename mut_index_cont,typename mcont_t,
	     typename mcounts_t,typename mutation_removal_policy>
    inline typename mut_index_cont::const_iterator find_fixation_pol(const mut_index_cont & mc,
								     const mcont_t & mutations,
								     const mcounts_t & mcounts,
								     const uint_t twoN,
								     mutation_removal_policy && mp) noexcept
    {
      return std::find_if(mc.cbegin(),mc.cend(),[&mutations,&mcounts,&twoN,&mp](const std::size_t & i) noexcept
			  {
			    return mcounts[i]==twoN && mp(mutations[i]);
			  });
    }

    //! Wrapper for erase/remove idiom with a custom policy
    template<typename mut_index_cont,
	     typename mcont_t,
	     typename mcounts_t,
	     typename mutation_removal_policy>
    inline void gamete_cleaner_erase_remove_idiom_wrapper_pol( mut_index_cont & mc,
							       const mcont_t & mutations,
							       const mcounts_t & mcounts,
							       const typename mut_index_cont::value_type & first_fixation,
							       const uint_t twoN,
							       mutation_removal_policy && mp) noexcept
    {
      /*
	The first call to std::find relies on a de-referencing of the return value of 
	find_fixation, passed to this function as first_fixation.  Because mc is sorted according to mutation position, 
	first_fixation is the fixation with the smallest position.  Using std::find like this avoids some calls to the lambda
	expression, which experiences cache-misses because mcounts is NOT sorted according to mutation position.
      */
      mc.erase(std::remove_if(std::find(mc.begin(),mc.end(),first_fixation),
			      mc.end(),[&mcounts,&mutations,&twoN,&mp](const typename mut_index_cont::value_type & i) noexcept
			      {
				return mcounts[i]==twoN && mp(mutations[i]);
			      }),mc.end());
    }

    /*! \brief Handles removal of indexes to mutations from gametes after sampling
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

    /*
      Now, we have the two overloads of gamete_cleaner.

      Future versions should consider merging into one policy-based function
      called by two different wrappers.  For now, they are kept separate.
    */
    
    /*! \brief Handles removal of indexes to mutations from gametes after sampling
      Intended use is when std::is_same< mutation_removal_policy, KTfwd::true_type >::type is true.
      Called by KTfwd::sample_diploid
    */
    template<typename gcont_t,typename mcont_t,typename mutation_removal_policy>
    inline typename std::enable_if< std::is_same<mutation_removal_policy,std::true_type>::value >::type
    gamete_cleaner(gcont_t & gametes, const mcont_t &, const std::vector<uint_t> & mcounts,
		   const uint_t twoN, const mutation_removal_policy &)
    {
      auto extant_gamete = next_extant_gamete(gametes.begin(),gametes.end());
      if(extant_gamete==gametes.end()) return;
      const auto fixation_n = find_fixation(extant_gamete->mutations,mcounts,twoN);
      bool neutral_fixations_exist = (fixation_n!=extant_gamete->mutations.cend());
      const auto fixation_s = find_fixation(extant_gamete->smutations,mcounts,twoN);
      bool selected_fixations_exist = (fixation_s!=extant_gamete->smutations.cend());
      if(!neutral_fixations_exist && !selected_fixations_exist) return;

      const auto gend = gametes.end();
      //Assign values to avoid tons of de-referencing later
      const auto fixation_n_value = (fixation_n==extant_gamete->mutations.cend()) ? typename decltype(fixation_n)::value_type() : *fixation_n;
      const auto fixation_s_value = (fixation_s==extant_gamete->smutations.cend()) ? typename decltype(fixation_s)::value_type() : *fixation_s;
      while(extant_gamete < gend)
	{

	  if(neutral_fixations_exist)
	    {
	      gamete_cleaner_erase_remove_idiom_wrapper(extant_gamete->mutations,mcounts,fixation_n_value,twoN);
	    }
	  if(selected_fixations_exist)
	    {
	      gamete_cleaner_erase_remove_idiom_wrapper(extant_gamete->smutations,mcounts,fixation_s_value,twoN);
	    }
	  extant_gamete = next_extant_gamete(extant_gamete+1,gend);
	}
    }

    /*! \brief Handles removal of indexes to mutations from gametes after sampling
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
      auto extant_gamete = next_extant_gamete(gametes.begin(),gametes.end());
      if(extant_gamete==gametes.end()) return;
      const auto fixation_n = find_fixation_pol(extant_gamete->mutations,mutations,mcounts,twoN,std::forward<decltype(mp)>(mp));
      bool neutral_fixations_exist = (fixation_n!=extant_gamete->mutations.cend());
      const auto fixation_s = find_fixation_pol(extant_gamete->smutations,mutations,mcounts,twoN,std::forward<decltype(mp)>(mp));
      bool selected_fixations_exist = (fixation_s!=extant_gamete->smutations.cend());
      if(!neutral_fixations_exist && !selected_fixations_exist) return;

      const auto gend = gametes.end();
      //Assign values to avoid tons of de-referencing later
      const auto fixation_n_value = (fixation_n==extant_gamete->mutations.cend()) ? typename decltype(fixation_n)::value_type() : *fixation_n;
      const auto fixation_s_value = (fixation_s==extant_gamete->smutations.cend()) ? typename decltype(fixation_s)::value_type() : *fixation_s;
      while(extant_gamete < gend)
	{
	  if(neutral_fixations_exist)
	    {
	      gamete_cleaner_erase_remove_idiom_wrapper_pol(extant_gamete->mutations,mutations,
							    mcounts,fixation_n_value,
							    twoN,std::forward<decltype(mp)>(mp));
	    }
	  if(selected_fixations_exist)
	    {
	      gamete_cleaner_erase_remove_idiom_wrapper_pol(extant_gamete->smutations,mutations,
							    mcounts,fixation_s_value,
							    twoN,std::forward<decltype(mp)>(mp));
	    }
	  extant_gamete = next_extant_gamete(extant_gamete+1,gend);
	}
    }
  }
}

#endif
