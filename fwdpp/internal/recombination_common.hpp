#ifndef __FDWPP_INTERNAL_RECOMBINATION_COMMON_HPP__
#define __FDWPP_INTERNAL_RECOMBINATION_COMMON_HPP__

#include <fwdpp/internal/rec_gamete_updater.hpp>
#include <cassert>
namespace KTfwd {

  namespace fwdpp_internal {

    template<typename double_vec_type,
	     typename gamete_type,
	     typename gamete_cont_iterator >
    void recombine_gametes( const double_vec_type & pos,
			    gamete_cont_iterator & ibeg,
			    gamete_cont_iterator & jbeg,
			    gamete_type & new_gamete1,
			    gamete_type & new_gamete2 )
    {
      assert( std::is_sorted(pos.cbegin(),pos.cend()) );
      short SWITCH = 0;

      auto itr = ibeg->mutations.cbegin(),
	jtr = jbeg->mutations.cbegin(),
	itr_s = ibeg->smutations.cbegin(),
	jtr_s = jbeg->smutations.cbegin(),
	itr_e = ibeg->mutations.cend(),
	itr_s_e = ibeg->smutations.cend(),
	jtr_e = jbeg->mutations.cend(),
	jtr_s_e = jbeg->smutations.cend();
#ifndef NDEBUG
      decltype(ibeg->mutations.size()) nm1=ibeg->mutations.size()+ibeg->smutations.size(),
	nm2=jbeg->mutations.size()+jbeg->smutations.size();
#endif
      for( const auto dummy : pos )
	{
	  itr = fwdpp_internal::rec_gam_updater(itr,itr_e,
						new_gamete2.mutations,new_gamete1.mutations,SWITCH,dummy);
	  itr_s = fwdpp_internal::rec_gam_updater(itr_s,itr_s_e,
						  new_gamete2.smutations,new_gamete1.smutations,SWITCH,dummy);
	  jtr = fwdpp_internal::rec_gam_updater(jtr,jtr_e,
						new_gamete1.mutations,new_gamete2.mutations,SWITCH,dummy);
	  jtr_s = fwdpp_internal::rec_gam_updater(jtr_s,jtr_s_e,
						  new_gamete1.smutations,new_gamete2.smutations,SWITCH,dummy);
	  SWITCH=!SWITCH;
	}
#ifndef NDEBUG
      decltype(new_gamete1.mutations.size()) __nm1 = new_gamete1.mutations.size()+new_gamete1.smutations.size(),
	__nm2 = new_gamete2.mutations.size()+new_gamete2.smutations.size();
      assert(__nm1+__nm2 == nm1+nm2);
#endif

      //Through fwdpp 0.2.4, we did a sort here, but it is not necessary
      /*
	std::sort(new_gamete1.smutations.begin(),new_gamete1.smutations.end(),
	[](typename gamete_type::mutation_list_type_iterator lhs,
	typename gamete_type::mutation_list_type_iterator rhs) { return lhs->pos < rhs->pos; });
	std::sort(new_gamete1.smutations.begin(),new_gamete1.smutations.end(),
	[](typename gamete_type::mutation_list_type_iterator lhs,
	typename gamete_type::mutation_list_type_iterator rhs) { return lhs->pos < rhs->pos; });
	std::sort(new_gamete2.mutations.begin(),new_gamete2.mutations.end(),
	[](typename gamete_type::mutation_list_type_iterator lhs,
	typename gamete_type::mutation_list_type_iterator rhs) { return lhs->pos < rhs->pos; });
	std::sort(new_gamete2.smutations.begin(),new_gamete2.smutations.end(),
	[](typename gamete_type::mutation_list_type_iterator lhs,
	typename gamete_type::mutation_list_type_iterator rhs) { return lhs->pos < rhs->pos; });
      */
#ifndef NDEBUG
      using mlist_itr = typename gamete_type::mutation_list_type_iterator;
      auto am_I_sorted = [](mlist_itr lhs,mlist_itr rhs){return lhs->pos < rhs->pos;};
      assert( std::is_sorted(new_gamete1.mutations.begin(),new_gamete1.mutations.end(),std::cref(am_I_sorted)) );
      assert( std::is_sorted(new_gamete1.smutations.begin(),new_gamete1.smutations.end(),std::cref(am_I_sorted)) );
      assert( std::is_sorted(new_gamete2.mutations.begin(),new_gamete2.mutations.end(),std::cref(am_I_sorted)) );
      assert( std::is_sorted(new_gamete2.smutations.begin(),new_gamete2.smutations.end(),std::cref(am_I_sorted)) );
#endif
    }
    
  }
}
#endif
