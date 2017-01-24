#ifndef __FWDPP_SUGAR_GSL_DELETER_HPP__
#define __FWDPP_SUGAR_GSL_DELETER_HPP__

#include <gsl/gsl_rng.h>
#include <memory>

namespace KTfwd
{
    namespace sugar
    {
        /*!
          \brief Function object to delete a gsl_rng *
          \ingroup sugar
         */
        struct gsl_rng_deleter
        {
            //! calls gsl_rng_free on r
            void
            operator()(gsl_rng *r) noexcept
            {
                gsl_rng_free(r);
            }
        };

        /*!
          \brief Smart pointer wrapper to gsl_rng *
          \ingroup sugar
        */
        using gsl_rng_ptr_t = std::unique_ptr<gsl_rng, gsl_rng_deleter>;
    }
}

#endif
