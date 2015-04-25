#ifndef __FWDPP_SUGAR_MUTATION_POPGENMUT_HPP__
#define __FWDPP_SUGAR_MUTATION_POPGENMUT_HPP__

#include <fwdpp/forward_types.hpp>

namespace KTfwd
{
  /*!
    The "standard" mutation type for population genetic simulation.
    A mutation has its own selection and dominance coefficients.
   */
  struct popgenmut : public mutation_base
  {
    //! The generation when the mutation arose
    unsigned g;
    //! Selection and dominance coefficients
    double s,h;
    popgenmut(const double & __pos, const double & __s, const double & __h,
	      const unsigned & __g,const unsigned & __n)
      : mutation_base(__pos,__n,(__s==0.) ? true : false),g(__g),s(__s),h(__h)
    {	
    }
  };
}
#endif
