/*! \file fwdpp/sugar/gsl/tags.hpp
 */
#ifndef __FWDPP_SUGAR_GSL_TAGS_HPP__
#define __FWDPP_SUGAR_GSL_TAGS_HPP__

namespace KTfwd
{
    namespace sugar
    {
        /*!
          \brief gsl_rng * types supported by fwdpp's sugar layer
          \ingroup sugar
         */
        enum class GSL_RNG_TYPE
        {
            MT19937,
            TAUS2
        };
        /*!
          \brief Dispatch tag for gsl_rng * types
          \ingroup sugar
        */
        template <GSL_RNG_TYPE> struct GSL_RNG_TYPE_TAG
        {
        };
    }
}

#endif
