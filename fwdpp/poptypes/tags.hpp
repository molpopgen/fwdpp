/*! \file fwdpp/sugar/poptypes/tags.hpp
 */
#ifndef __FWDPP_SUGAR_POPTYPES_TAGS_HPP__
#define __FWDPP_SUGAR_POPTYPES_TAGS_HPP__

namespace fwdpp
{
    namespace sugar
    {
        //! Types of populations supported by sugar layer \ingroup sugar
        enum class FWDPP_SUGAR_POPTYPE
        {
            SINGLELOC,
            MULTILOC
        };
        //! Dispatch tag template for population types supported by sugar layer
        //! \ingroup sugar
        template <FWDPP_SUGAR_POPTYPE> struct FWDPP_SUGAR_POPTAG
        {
        };
        //! Single-locus simulations \ingroup sugar
        using SINGLELOC_TAG = FWDPP_SUGAR_POPTAG<FWDPP_SUGAR_POPTYPE::SINGLELOC>;
        //! Multi-locus simulations \ingroup sugar
        using MULTILOC_TAG
            = FWDPP_SUGAR_POPTAG<FWDPP_SUGAR_POPTYPE::MULTILOC>;
    }
}

#endif
