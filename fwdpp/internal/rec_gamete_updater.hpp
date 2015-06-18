#ifndef __FWDPP_INTERNAL_REC_GAMETE_UPDATER_HPP__
#define __FWDPP_INTERNAL_REC_GAMETE_UPDATER_HPP__

#include <algorithm>
#include <iterator>
#include <functional>

namespace KTfwd
{
  namespace fwdpp_internal
  {
    template< typename itr_type,
	      typename cont_type >
    itr_type rec_gam_updater( itr_type __first, itr_type __last,
			      cont_type & muts,
			      const double & val )
    {
      //O(log_2) comparisons of double plus at most __last - __first copies
      itr_type __ub = std::lower_bound(__first,__last,
				       std::cref(val),
				       [](const typename itr_type::value_type & __mut,const double & __val) {
					 return __mut->pos < __val;
				       });
      /*
	NOTE: the use of insert here
	instead of std::copy(__first,__ub,std::back_inserter(muts));
	Reduced peak RAM use on GCC 4.9.2/Ubuntu Linux.
      */
      muts.insert(muts.end(),__first,__ub);
      return __ub;
    }
  }
}

#endif
