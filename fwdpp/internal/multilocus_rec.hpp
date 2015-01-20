#ifndef __FWDPP_INTERNAL_MULTILOCUS_REC_HPP__
#define __FWDPP_INTERNAL_MULTILOCUS_REC_HPP__

/*
  The mechanics of crossing over for a multilocus
  simulation
 */

#include <gsl/gsl_rng.h>

namespace KTfwd {
  namespace fwdpp_internal {
    template<typename within_loc_rec_policy,
	     typename between_loc_rec_policy,
	     typename gamete_itr_t>
    gamete_itr_t multilocus_rec(gsl_rng * r,
				const within_loc_rec_policy & rec,
				const between_loc_rec_policy & bw,
				const double * r_between_loci,
				const unsigned & i,
				gamete_itr_t & parental_gamete_1,
				gamete_itr_t & parental_gamete_2,
				bool & g1, bool & LO )
    {
      unsigned temp = rec( parental_gamete_1,parental_gamete_2 );
      if ( i > 0 )
	{
	  unsigned nrbw = bw(r,r_between_loci[i-1]);
	  bool obw = (nrbw%2!=0) ? true : false;
	  g1 = (LO) ? !g1 : g1;
	  g1 = (obw) ? !g1 : g1;
	}
      LO = (temp % 2 != 0.) ? true : false;
      return (g1) ? parental_gamete_1 : parental_gamete_2;
    }
  }
}

#endif
