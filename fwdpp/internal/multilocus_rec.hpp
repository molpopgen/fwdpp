#ifndef __FWDPP_INTERNAL_MULTILOCUS_REC_HPP__
#define __FWDPP_INTERNAL_MULTILOCUS_REC_HPP__

/*
  The mechanics of crossing over for a multilocus
  simulation
 */

#include <gsl/gsl_rng.h>
#include <iostream>
namespace KTfwd {
  namespace fwdpp_internal {
    /*!
      The mechanics of recombination for multilocus simulations
      \param r A GSL random number generator
      \param rec A policy that affects recombination withing a locus
      \param bw A policy that returns the number of recombination events between loci.  The return value must be unsigned int, and what really matters is if that value is odd or even.
      \param r_between_loci A const array of doubles corresponding to the recombination rates between loci.  For k loci, this vector must be k-1 doubles long.  Further, the i-th value must correspond to the recombination rate between loci i-1 and i.
      \param parental_gamete_1 An iterator to the first parental gamete
      \param parental_gamete_1 An iterator to the second parental gamete
      \param g1 If true, then parental_gamete_1 is inherited by the descendant.
      \paragm LO Should be true if the last recombination even within loci resulted in an odd number of crossovers.
     */
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
      std::cerr << "mloc_internal: " << parental_gamete_1->mutations.size() << ' ' << parental_gamete_2->mutations.size() << ' ' << temp << '\n';
      if ( i > 0 )
	{
	  std::cerr << "HERE!\n";
	  unsigned nrbw = bw(r,r_between_loci[i-1]);
	  bool obw = (nrbw%2!=0) ? true : false;
	  g1 = (LO) ? !g1 : g1;
	  g1 = (obw) ? !g1 : g1;
	}
      LO = (temp % 2 != 0.) ? true : false;
      std::cerr << "g1 = " << g1 << ' ' << LO << '\n';
      return (g1) ? parental_gamete_1 : parental_gamete_2;
    }
  }
}

#endif
