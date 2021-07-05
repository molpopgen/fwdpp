#pragma once

#include <type_traits>

namespace fwdpp
{
    namespace ts
    {
        namespace types
        {
            template <typename SignedInteger>
            inline constexpr SignedInteger
            generate_null_id()
            {
                static_assert(std::is_integral<SignedInteger>::value,
                              "T must be an integral type");
                static_assert(std::is_signed<SignedInteger>::value,
                              "T must be a signed type");
                return SignedInteger{-1};
            }
        }
    }
}
