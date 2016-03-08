#ifndef __FDWPP_INTERNAL_RECOMBINATION_COMMON_HPP__
#define __FDWPP_INTERNAL_RECOMBINATION_COMMON_HPP__

#include <fwdpp/internal/rec_gamete_updater.hpp>
#include <cassert>
namespace KTfwd {

  namespace fwdpp_internal {

    template<typename double_vec_type,
	     typename gcont_t,
	     typename mcont_t>
    void recombine_gametes( const double_vec_type & pos,
			    const std::size_t ibeg,
			    const std::size_t jbeg,
			    gcont_t & gametes,
			    const mcont_t & mutations,
			    typename gcont_t::value_type::mutation_container & neutral,
			    typename gcont_t::value_type::mutation_container & selected)
    {
      assert( std::is_sorted(pos.cbegin(),pos.cend()) );

      auto itr = gametes[ibeg].mutations.cbegin();
      auto jtr = gametes[jbeg].mutations.cbegin();
      auto itr_s = gametes[ibeg].smutations.cbegin();
      auto jtr_s = gametes[jbeg].smutations.cbegin();
      auto itr_e = gametes[ibeg].mutations.cend();
      auto itr_s_e = gametes[ibeg].smutations.cend();
      auto jtr_e = gametes[jbeg].mutations.cend();
      auto jtr_s_e = gametes[jbeg].smutations.cend();

      auto pend = (pos.size()%2==0) ? pos.cend() : pos.cend()-1;
      auto i = pos.cbegin();
      const double * p = pos.data();
      for( ; i < pend ; i+=2,p+=2 )
	{
	  itr = fwdpp_internal::rec_gam_updater(itr,itr_e,mutations,
						neutral,*p);
	  itr_s = fwdpp_internal::rec_gam_updater(itr_s,itr_s_e,mutations,
						  selected,*p);
	  jtr = fwdpp_internal::rec_update_itr(jtr,jtr_e,mutations,*p);
	  jtr_s = fwdpp_internal::rec_update_itr(jtr_s,jtr_s_e,mutations,*p);

	  jtr = fwdpp_internal::rec_gam_updater(jtr,jtr_e,mutations,
						neutral,*(p+1));
	  jtr_s = fwdpp_internal::rec_gam_updater(jtr_s,jtr_s_e,mutations,
						  selected,*(p+1));
	  itr = fwdpp_internal::rec_update_itr(itr,itr_e,mutations,*(p+1));
	  itr_s = fwdpp_internal::rec_update_itr(itr_s,itr_s_e,mutations,*(p+1));
	}
      assert(std::distance(i,pos.cend())<=1);
      for( ; i<pos.cend();++i)
	{
	  itr = fwdpp_internal::rec_gam_updater(itr,itr_e,mutations,
						neutral,*p);
	  itr_s = fwdpp_internal::rec_gam_updater(itr_s,itr_s_e,mutations,
						  selected,*p);
	}
      assert(gamete_is_sorted_n(gametes[ibeg],mutations));
      assert(gamete_is_sorted_s(gametes[ibeg],mutations));
    }    
  }
}
#endif
