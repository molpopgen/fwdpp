#ifndef FWDPP_META_ALWAYS_FALSE_HPP__
#define FWDPP_META_ALWAYS_FALSE_HPP__

#include <type_traits>

namespace fwdpp
{
    namespace meta
    {
        template <typename T> struct always_false : std::false_type
        {
        };
    }
}

#endif
