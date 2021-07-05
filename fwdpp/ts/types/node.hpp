#pragma once

#include <tuple>
#include <type_traits>
#include "generate_null_id.hpp"

namespace fwdpp
{
    namespace ts
    {
        namespace types
        {
            template <typename SignedInteger> struct node
            /// A node in a tree sequence
            ///  \version 0.7.0 Added to fwdpp
            ///  \version 0.8.0 Rename ::population to ::deme
            {
                static_assert(std::is_integral<SignedInteger>::value,
                              "SignedInteger must be an integral type");
                static_assert(std::is_signed<SignedInteger>::value,
                              "SignedInteger must be a signed type");
                using id_type = SignedInteger;
                /// Null ID value
                static constexpr SignedInteger null = generate_null_id<id_type>();
                /// Location of the node.
                /// Used for models of discrete population structure
                id_type deme;
                /// Birth time of the node.
                double time;
            };

#if __cplusplus < 201703L
            template <typename SignedInteger>
            constexpr SignedInteger node<SignedInteger>::null;
#endif

            template <typename SignedInteger>
            inline bool
            operator==(const node<SignedInteger>& self,
                       const node<SignedInteger>& other)
            {
                return std::tie(self.time, self.deme)
                       == std::tie(other.time, other.deme);
            }
        }
    }
}
