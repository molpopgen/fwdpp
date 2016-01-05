#ifndef __FWDPP_GENETICS101_HPP__
#define __FWDPP_GENETICS101_HPP__

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

namespace KTfwd
{
  struct genetics101
  /*! Genetics 101: simple model of recombination.  r is the probability that the two gametes recombine
   */
  {
    using result_type = unsigned;
    template<typename glookup_t,
	     typename queue_t,
	     typename gcont_t,
	     typename mcont_t,
	     typename rec_pos_generator>
    std::size_t operator()( const size_t g1,
			    const size_t g2,
			    glookup_t & gamete_lookup,
			    queue_t & gamete_recycling_bin,
			    std::vector<std::size_t> & neutral,
			    std::vector<std::size_t> & selected,
			    gcont_t & gametes,
			    const mcont_t & mutations,
			    const double & littler,
			    gsl_rng * r,
			    const rec_pos_generator & rp) const
    {
      return recombine_gametes(r,littler,gametes,mutations,g1,g2,gamete_lookup,gamete_recycling_bin,neutral,selected,rp);
    }
  };
}

#endif //__FWDPP_GENETICS101_HPP__
