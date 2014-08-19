#ifndef __FWDPP_INTERNAL_REC_GAMETE_UPDATER_HPP__
#define __FWDPP_INTERNAL_REC_GAMETE_UPDATER_HPP__

namespace KTfwd
{
  namespace fwdpp_internal
  {
    struct rec_gamete_updater
    {
      template<typename itr_type,
	       typename cont_type>
      inline bool operator()( itr_type & i, cont_type * m1, cont_type * m2,
			      const short & SWITCH, const double & val ) const
      {
	if( i->pos < val )
	  {
	    if( SWITCH )
	      {
		m1->emplace_back(i);
	      }
	    else
	      {
		m2->emplace_back(i);
	      }
	    return true;
	  }
	return false;
      }
    };
  }
}

#endif
