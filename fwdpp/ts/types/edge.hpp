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
            template <typename SignedInteger> struct edge
            /*! An edge in a tree sequence
             * @version 0.7.0 Added to fwdpp
             *
             * Edges define a transmission event
             * of the genomic interval [left,right)
             * from parent to child.
             */
            {
                static_assert(std::is_integral<SignedInteger>::value,
                              "SignedInteger must be an integral type");
                static_assert(std::is_signed<SignedInteger>::value,
                              "SignedInteger must be a signed type");
                using id_type = SignedInteger;
                /// Null ID value
                static constexpr SignedInteger null = generate_null_id<id_type>();
                /// Left (inclusive) edge of genomic segment
                double left;
                /// Right (exclusive) edge of genomic segment
                double right;
                /// Parent ID
                id_type parent;
                /// Child ID
                id_type child;
            };

#if __cplusplus < 201703L
            template <typename SignedInteger>
            constexpr SignedInteger edge<SignedInteger>::null;
#endif

            template <typename SignedInteger>
            inline bool
            operator==(const edge<SignedInteger>& self,
                       const edge<SignedInteger>& other)
            {
                return std::tie(self.parent, self.child, self.left, self.right)
                       == std::tie(other.parent, other.child, other.left, other.right);
            }
        }
    }
}
