#ifndef __FWDPP_SUGAR_GSL_DELETER_HPP__
#define __FWDPP_SUGAR_GSL_DELETER_HPP__

#include <gsl/gsl_rng.h>
#include <memory>

namespace fwdpp
{
    namespace gsl
    {
        /*!
          \brief Smart pointer wrapper to gsl_rng *
          \ingroup sugar
        */
        using gsl_rng_ptr_t = std::unique_ptr<gsl_rng, void (*)(gsl_rng *)>;
    }
}

#endif
