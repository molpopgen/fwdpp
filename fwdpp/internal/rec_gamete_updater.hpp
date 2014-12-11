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
			      cont_type & m1, cont_type & m2,
			      const short & SWITCH, const double & val )
    {
      typename itr_type::difference_type __dist = std::distance(__first,__last);
      if ( __dist <= 250 )
	{
	  //At most __last - _first comparisons of double + equal number if emplace_backs
	  while(__first != __last && (*__first)->pos < val)
	    {
	      if( SWITCH )
		{
		  m1.emplace_back(*__first);
		}
	      else
		{
		  m2.emplace_back(*__first);
		}
	      ++__first;
	    }
	  return __first;
	}
      //O(log_2) comparisons of double plus at most __last - __first copies
      itr_type __ub = std::upper_bound(__first,__last,
				       std::cref(val),
				       [](const double __val,const typename itr_type::value_type __mut) {
					 return __val < __mut->pos;
				       });
      if (SWITCH)
	std::copy(__first,__ub,std::back_inserter(m1));
      else
	std::copy(__first,__ub,std::back_inserter(m2));
      return __ub;
    }
  }
}

#endif
