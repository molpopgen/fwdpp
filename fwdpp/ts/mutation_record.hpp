#ifndef FWDPP_TS_MUTATION_RECORD_HPP
#define FWDPP_TS_MUTATION_RECORD_HPP

#include <tuple>
#include <cstdint>

namespace fwdpp
{
    namespace ts
    {
        struct mutation_record
        /// Tracks mutations on tree sequences.
        ///  \version 0.7.0 Added to fwdpp
        {
            /// The node to which the mutation is
            /// currently simplified
            std::int32_t node;
            /// The index of the mutation in the
            /// population's mutation container
            std::size_t key;
            /// Row in the site table.
            std::size_t site_id;
            /// Character state of the mutation
            std::int8_t derived_state; // TODO: should this be a template type?
            /// True if mutation affects fitness, otherwise false.
            bool neutral;
        };

        inline bool
        operator==(const mutation_record& a, const mutation_record& b)
        {
            if (a.site_id == b.site_id)
                {
                    return std::tie(a.node, a.key, a.derived_state, a.neutral)
                           == std::tie(b.node, b.key, b.derived_state,
                                       b.neutral);
                }
            return false;
        }
    } // namespace ts
} // namespace fwdpp

#endif
