/*
  This file provides lightweight wrappers
  around calls to GSL.

  The intendend use is to implement 
  standard models of distributions on selection/
  dominance effects of mutations.

  Rather, the intended use is that they provide
  the implementation details of "shmodel", defined below
*/

#ifndef __FWDPP_EXTENSIONS_CALLBACKS_HPP__
#define __FWDPP_EXTENSIONS_CALLBACKS_HPP__

#include <gsl/gsl_randist.h>
#include <functional>
#include <exception>
#include <cmath>

namespace KTfwd {

  namespace extensions {

    struct shmodel
    {
      std::function<double(gsl_rng*)> s,h;
    };
  
    struct constant
    {
      double x;
      constant(const double & __x) : x(__x)
      {
	if(!std::isfinite(x)) {
	  throw std::runtime_error("value must be finite");
	}
      }
      inline double operator()(gsl_rng *) const
      {
	return x;
      }
    };

    struct exponential
    {
      double mean;
      exponential(const double & m) : mean(m)
      {
	if(!std::isfinite(mean))
	  {
	    throw std::runtime_error("mean must be finite");
	  }
	if(mean==0.)
	  {
	    throw std::runtime_error("mean must not equal 0");
	  }
      }
      inline double operator()(gsl_rng * r) const
      {
	return gsl_ran_exponential(r,mean);
      }
    };

    struct uniform
    {
      double mn,mx;
      uniform(const double & __mn,
	      const double & __mx) : mn(__mn),mx(__mx)
      {
	if(!std::isfinite(mn) || !std::isfinite(mx))
	  {
	    throw std::runtime_error("min and max of range must both be finite");
	  }
	if(mn>mx)
	  {
	    throw std::runtime_error("min must be <= max");
	  }
      }
      inline double operator()(gsl_rng * r) const
      {
	return gsl_ran_flat(r,mn,mx);
      }
    };

    struct beta
    {
      double a,b,factor;
      beta(const double & __a,
	   const double & __b,
	   const double & __f) : a(__a),b(__b),factor(__f)
      {
	if(!std::isfinite(a) || a <= 0.)
	  {
	    throw std::runtime_error("a must be > 0.");
	  }
	if(!std::isfinite(b) || b <= 0.)
	  {
	    throw std::runtime_error("b must be > 0.");
	  }
	if(!std::isfinite(factor) || !(factor>0.))
	  {
	    throw std::runtime_error("scaling factor must be finite and > 0");
	  }
      }
      inline double operator()(gsl_rng * r) const
      {
	return factor*gsl_ran_beta(r,a,b);
      }
    };

    struct gaussian
    {
      double sd;
      gaussian(const double & __sd) : sd(__sd)
      {
	if(sd == 0.) throw std::runtime_error("sd must not equal 0");
	if(!std::isfinite(sd)) throw std::runtime_error("sd must be finite");
      }
      inline double operator()(const gsl_rng * r) const
      {
	return gsl_ran_gaussian(r,sd);
      }
    };

    struct gamma
    {
      double mean,shape;
      gamma(const double & __m,
	    const double & __s ) : mean(__m),shape(__s)
      {
	if(!std::isfinite(mean))
	  {
	    throw std::runtime_error("mean and sign must both be finite");
	  }
	if(shape <= 0)
	  {
	    throw std::runtime_error("shape must be >= 0");
	  }
      }
      inline double operator()(const gsl_rng * r) const
      {
	return gsl_ran_gamma(r,shape,mean/shape);
      }
    };
  }
}
#endif
