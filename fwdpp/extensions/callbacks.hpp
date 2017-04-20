/*!
  \file callbacks.hpp

  This file provides lightweight wrappers
  around calls to GSL.

  The intended use is to implement
  standard models of distributions on selection/
  dominance effects of mutations.

  Rather, the intended use is that they provide
  the implementation details of "shmodel", defined below
*/

#ifndef __FWDPP_EXTENSIONS_CALLBACKS_HPP__
#define __FWDPP_EXTENSIONS_CALLBACKS_HPP__

#include <gsl/gsl_randist.h>
#include <functional>
#include <stdexcept>
#include <cmath>

namespace KTfwd
{

    /*!
      Useful types for implementing fwdpp-based simulations in
      enviroments like R, Python, etc.

      Unlike the rest of fwdpp, functions in this namespace
      are allowed to throw exceptions, and it is up to the programmer
      to catch them and handle them appropriately.  Both Rcpp and
      Cython/boost.python

      Examples of using this namespace are:

      1. http://github.com/molpopgen/fwdpy (This project is deprecated
      and no longer tracks changes to the extensions API.)
      2. http://github.com/molpopgen/fwdpy11
     */
    namespace extensions
    {

        struct shmodel
        /*!
          Callback wrapper.  Used to model
          distributions on effect sizes/selection coefficients (s)
          and dominance (h)
        */
        {
            std::function<double(const gsl_rng *)> s, h;
            //! Default constructor useful in extension situations that don't
            //! understand std::function
            shmodel() = default;
            //! More efficient constructor for c++11-aware situations
            shmodel(std::function<double(const gsl_rng *)> sfxn,
                    std::function<double(const gsl_rng *)> hfxn)
                : s(std::move(sfxn)), h(std::move(hfxn))
            {
            }
        };

        struct constant
        /*!
          Callback for fixed s and/or h
         */
        {
            const double x;
            constant(const double &__x) : x(__x)
            {
                if (!std::isfinite(x))
                    {
                        throw std::invalid_argument("value must be finite");
                    }
            }
            inline double
            operator()(const gsl_rng *) const
            {
                return x;
            }
        };

        struct exponential
        /*!
          Exponential s or h
         */
        {
            const double mean;
            exponential(const double &m) : mean(m)
            {
                if (!std::isfinite(mean))
                    {
                        throw std::invalid_argument("mean must be finite");
                    }
                if (mean == 0.)
                    {
                        throw std::invalid_argument("mean must not equal 0");
                    }
            }
            inline double
            operator()(const gsl_rng *r) const
            {
                return gsl_ran_exponential(r, mean);
            }
        };

        struct uniform
        /*!
          Uniform s or h
         */
        {
            const double mn, mx;
            uniform(const double &__mn, const double &__mx)
                : mn(__mn), mx(__mx)
            {
                if (!std::isfinite(mn) || !std::isfinite(mx))
                    {
                        throw std::invalid_argument(
                            "min and max of range must both be finite");
                    }
                if (mn > mx)
                    {
                        throw std::invalid_argument("min must be <= max");
                    }
            }
            inline double
            operator()(const gsl_rng *r) const
            {
                return gsl_ran_flat(r, mn, mx);
            }
        };

        struct beta
        /*!
          Beta-distributed s or h
        */
        {
            const double a, b, factor;
            beta(const double &__a, const double &__b, const double &__f = 1)
                : a(__a), b(__b), factor(__f)
            {
                if (!std::isfinite(a) || a <= 0.)
                    {
                        throw std::invalid_argument("a must be > 0.");
                    }
                if (!std::isfinite(b) || b <= 0.)
                    {
                        throw std::invalid_argument("b must be > 0.");
                    }
                if (!std::isfinite(factor) || !(factor > 0.))
                    {
                        throw std::invalid_argument(
                            "scaling factor must be finite and > 0");
                    }
            }
            inline double
            operator()(const gsl_rng *r) const
            {
                return factor * gsl_ran_beta(r, a, b);
            }
        };

        struct gaussian
        /*!
          Gaussian s or h
        */
        {
            const double sd;
            gaussian(const double &__sd) : sd(__sd)
            {
                if (!(sd > 0.))
                    throw std::invalid_argument("sd must be > 0");
                if (!std::isfinite(sd))
                    throw std::invalid_argument("sd must be finite");
            }
            inline double
            operator()(const gsl_rng *r) const
            {
                return gsl_ran_gaussian_ziggurat(r, sd);
            }
        };

        struct gamma
        /*!
          Gamma distributed s or h
        */
        {
            const double mean, shape;
            gamma(const double &__m, const double &__s) : mean(__m), shape(__s)
            {
                if (!std::isfinite(mean) || !std::isfinite(shape))
                    {
                        throw std::invalid_argument(
                            "mean and shape must both be finite");
                    }
                if (!(shape > 0))
                    {
                        throw std::invalid_argument("shape must be > 0");
                    }
            }
            inline double
            operator()(const gsl_rng *r) const
            {
                return gsl_ran_gamma(r, shape, mean / shape);
            }
        };
    }
}
#endif
