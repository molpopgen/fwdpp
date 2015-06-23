#ifndef __FWDPP_INTERNAL_GAMETE_CLEANER_HPP__
#define __FWDPP_INTERNAL_GAMETE_CLEANER_HPP__

#include <type_traits>

namespace KTfwd {
  namespace fwdpp_internal {
    /*! \brief Handles removal of pointers to mutations from gametes after sampling
      Intended use is when std::is_same< mutation_removal_policy, KTfwd::remove_nothing >::type is true.
      Called by KTfwd::sample_diploid via dispatch.
    */
    template<typename gamete_list_type, typename mutation_removal_policy>
    inline void gamete_cleaner(gamete_list_type *, const mutation_removal_policy &, std::true_type) 
    {
      return;
    }

    /*! \brief Handles removal of pointers to mutations from gametes after sampling
      Intended use is when std::is_same< mutation_removal_policy, KTfwd::remove_nothing >::type is false.
      Called by KTfwd::sample_diploid via dispatch.
    */
    template<typename gamete_list_type, typename mutation_removal_policy>
    inline void gamete_cleaner(gamete_list_type * gametes, const mutation_removal_policy & mp, std::false_type) 
    {
      std::for_each( gametes->begin(),
		     gametes->end(),
		     [&mp]( typename gamete_list_type::value_type & __g ) {
		       __g.mutations.erase(std::remove_if(__g.mutations.begin(),__g.mutations.end(),std::cref(mp)),__g.mutations.end());
		       __g.smutations.erase(std::remove_if(__g.smutations.begin(),__g.smutations.end(),std::cref(mp)),__g.smutations.end());
		     });
    }
  }
}

#endif
