/*! \file fwdpp/poptypes/tags.hpp
 */
#ifndef __FWDPP_POPTYPES_TAGS_HPP__
#define __FWDPP_POPTYPES_TAGS_HPP__

namespace fwdpp
{
    namespace poptypes
    {
        //! Types of populations supported.
        enum class FWDPP_POPTYPE
        {
            DIPLOID
        };

        //! Dispatch tag template for population types supported by sugar layer
        template <FWDPP_POPTYPE> struct FWDPP_POPTAG
        {
        };
        //! Diploid simulations
        using DIPLOID_TAG = FWDPP_POPTAG<FWDPP_POPTYPE::DIPLOID>;
    }
}

#endif
