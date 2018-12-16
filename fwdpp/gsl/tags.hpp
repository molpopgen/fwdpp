/*! \file fwdpp/gsl/tags.hpp
 */
#ifndef __FWDPP_GSL_TAGS_HPP__
#define __FWDPP_GSL_TAGS_HPP__

namespace fwdpp
{
    namespace gsl
    {
        /*!
          \brief gsl_rng * types supported by fwdpp
         */
        enum class GSL_RNG_TYPE
        {
            MT19937,
            TAUS2
        };
        /*!
          \brief Dispatch tag for gsl_rng * types
        */
        template <GSL_RNG_TYPE> struct GSL_RNG_TYPE_TAG
        {
        };
    }
}

#endif
