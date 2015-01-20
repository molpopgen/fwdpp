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
      //unsigned temp = rec_policies[i]( p1c[i].first, p1c[i].second );
      unsigned temp = rec( parental_gamete_1,parental_gamete_2 );
      if ( i > 0 )
	{
	  unsigned nrbw = bw(r,r_between_loci[i-1]);
	  bool obw = (nrbw%2!=0) ? true : false;
	  g1 = (LO) ? !g1 : g1;
	  g1 = (obw) ? !g1 : g1;
	}
      //(ptr2cdip+i)->first = (p1g1) ? p1c[i].first : p1c[i].second;
      //offspring_gamete = (g1) ? parental_gamete_1 : parental_gamete_2;
      LO = (temp % 2 != 0.) ? true : false;
      return (g1) ? parental_gamete_1 : parental_gamete_2;
      // temp = rec_policies[i]( p2c[i].first, p2c[i].second );
      // if ( i > 0 )
      //   {
      //     unsigned nrbw = blrf(r,r_between_loci[i-1]);
      //     bool obw = (nrbw%2!=0) ? true : false;
      //     p2g1 = (LO2) ? !p2g1 : p2g1;
      //     p2g1 = (obw) ? !p2g1 : p2g1;
      //   }
      // (ptr2cdip+i)->second = (p2g1) ? p2c[i].first : p2c[i].second;
      // LO2 = (temp % 2 != 0.) ? true : false;
    }
  }
}

#endif
