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
        };

        bool
        operator==(const mutation_record& a, const mutation_record& b)
        {
            return std::tie(a.node, a.key) == std::tie(b.node, b.key);
        }
    } // namespace ts
} // namespace fwdpp

#endif
