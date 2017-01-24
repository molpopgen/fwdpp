#ifndef __FWDPP_TAGS_FITNESS_TAGS_HPP__
#define __FWDPP_TAGS_FITNESS_TAGS_HPP__

namespace KTfwd
{
    namespace tags
    {
        //! Dispatch tags for fitness models
        template <bool> struct diploid_type
        {
        };
        //! The default (first and second are passed on)
        using standard_diploid_t = diploid_type<false>;
        //! For more complex models that require all of the data in a custom
        //! diploid type
        using custom_diploid_t = diploid_type<true>;
    }
}

#endif
