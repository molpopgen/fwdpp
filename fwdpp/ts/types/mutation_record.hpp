#pragma once

#include <cstdint>
#include <tuple>
#include <type_traits>
#include "generate_null_id.hpp"

namespace fwdpp
{
    namespace ts
    {
        namespace types
        {
            template <typename SignedInteger> struct mutation_record
            /// Tracks mutations on tree sequences.
            ///  \version 0.7.0 Added to fwdpp
            {
                using id_type = SignedInteger;
                /// Null ID value
                static constexpr SignedInteger null = generate_null_id<id_type>();
                /// The node to which the mutation is
                /// currently simplified
                id_type node;
                /// The index of the mutation in the
                /// population's mutation container
                std::size_t key;
                /// Row in the site table.
                std::size_t site;
                /// Character state of the mutation
                std::int8_t derived_state; // TODO: should this be a template type?
                /// True if mutation affects fitness, otherwise false.
                bool neutral;
            };

#if __cplusplus < 201703L
            template <typename SignedInteger>
            constexpr SignedInteger mutation_record<SignedInteger>::null;
#endif

            template <typename SignedInteger>
            inline bool
            operator==(const mutation_record<SignedInteger>& self,
                       const mutation_record<SignedInteger>& other)
            {
                // Test site first b/c two mutations cannot be equal
                // if they aren't at the same site.
                return self.site == other.site
                       && std::tie(self.node, self.key, self.derived_state, self.neutral)
                              == std::tie(other.node, other.key, other.derived_state,
                                          other.neutral);
            }

        }
    }
}
