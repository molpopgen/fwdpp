#ifndef FWDPP_UTIL_ENUM_BITFLAGS_HPP
#define FWDPP_UTIL_ENUM_BITFLAGS_HPP

#include <type_traits>

namespace fwdpp
{
    namespace enum_bitflags
    {
        template <typename E>
        auto
        operator|(E a, E b)
        {
            using T = std::underlying_type_t<E>;
            static_assert(std::is_integral<T>::value,
                          "underlying type must be integral");
            static_assert(std::is_unsigned<T>::value,
                          "underlying type must be unsigned");
            return static_cast<T>(a) | static_cast<T>(b);
        }

        template <typename E>
        auto&
        operator|=(E& a, E b)
        {
            a = static_cast<E>(a | b);
            return a;
        }

        template <typename E>
        auto
        operator&(E a, E b)
        {
            using T = std::underlying_type_t<E>;
            static_assert(std::is_integral<T>::value,
                          "underlying type must be integral");
            static_assert(std::is_unsigned<T>::value,
                          "underlying type must be unsigned");
            return static_cast<T>(a) & static_cast<T>(b);
        }

        template <typename E>
        auto&
        operator&=(E& a, E b)
        {
            a = static_cast<E>(a & b);
            return a;
        }

        template <typename E>
        auto
        operator^(E a, E b)
        {
            using T = std::underlying_type_t<E>;
            static_assert(std::is_integral<T>::value,
                          "underlying type must be integral");
            static_assert(std::is_unsigned<T>::value,
                          "underlying type must be unsigned");
            return static_cast<T>(a) ^ static_cast<T>(b);
        }

        template <typename E>
        auto&
        operator^=(E& a, E b)
        {
            a = static_cast<E>(a ^ b);
            return a;
        }

        template <typename E>
        auto
        operator~(E a)
        {
            using T = std::underlying_type_t<E>;
            static_assert(std::is_integral<T>::value,
                          "underlying type must be integral");
            static_assert(std::is_unsigned<T>::value,
                          "underlying type must be unsigned");
            return ~static_cast<T>(a);
        }
    }
}

#endif
