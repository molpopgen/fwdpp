/*! \file fwdpp/sugar/poptypes/tags.hpp
 */
#ifndef __FWDPP_SUGAR_POPTYPES_TAGS_HPP__
#define __FWDPP_SUGAR_POPTYPES_TAGS_HPP__

namespace KTfwd
{
    namespace sugar
    {
        //! Types of populations supported by sugar layer \ingroup sugar
        enum class FWDPP_SUGAR_POPTYPE
        {
            SINGLE,
            META,
            MULTILOC
        };
        //! Dispatch tag template for population types supported by sugar layer
        //! \ingroup sugar
        template <FWDPP_SUGAR_POPTYPE> struct FWDPP_SUGAR_POPTAG
        {
        };
        //! Single-population simulations \ingroup sugar
        using SINGLEPOP_TAG = FWDPP_SUGAR_POPTAG<FWDPP_SUGAR_POPTYPE::SINGLE>;
        //! Meta-population simulations \ingroup sugar
        using METAPOP_TAG = FWDPP_SUGAR_POPTAG<FWDPP_SUGAR_POPTYPE::META>;
        //! Single-population, multi-locus simulations \ingroup sugar
        using MULTILOCPOP_TAG
            = FWDPP_SUGAR_POPTAG<FWDPP_SUGAR_POPTYPE::MULTILOC>;
    }
}

#endif
