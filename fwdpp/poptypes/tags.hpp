/*! \file fwdpp/poptypes/tags.hpp
 */
#ifndef __FWDPP_POPTYPES_TAGS_HPP__
#define __FWDPP_POPTYPES_TAGS_HPP__

namespace fwdpp
{
    namespace poptypes
    {
        //! Types of populations supported by sugar layer \ingroup sugar
        enum class FWDPP_POPTYPE
        {
            SINGLELOC,
            MULTILOC
        };
        //! Dispatch tag template for population types supported by sugar layer
        template <FWDPP_POPTYPE> struct FWDPP_POPTAG
        {
        };
        //! Single-locus simulations
        using SINGLELOC_TAG = FWDPP_POPTAG<FWDPP_POPTYPE::SINGLELOC>;
        //! Multi-locus simulations
        using MULTILOC_TAG
            = FWDPP_POPTAG<FWDPP_POPTYPE::MULTILOC>;
    }
}

#endif
