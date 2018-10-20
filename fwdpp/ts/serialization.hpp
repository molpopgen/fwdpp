#ifndef FWDPP_TS_SERIALIZATION_HPP
#define FWDPP_TS_SERIALIZATION_HPP

#include <cstdint>
#include <stdexcept>
#include <fwdpp/io/scalar_serialization.hpp>
#include "table_collection.hpp"

namespace fwdpp
{
    namespace ts
    {
        constexpr static const std::uint32_t TS_TABLES_VERSION = 1;

        template <typename ostreamtype, std::size_t version>
        struct serialize_node
        {
            inline void
            operator()(ostreamtype& o, const node& n) const
            {
                throw std::runtime_error("invalid serialization version");
            }
        };

        template <typename ostreamtype>
        struct serialize_node<ostreamtype, TS_TABLES_VERSION>
        {
            inline void
            operator()(ostreamtype& o, const node& n) const
            {
                io::scalar_writer sw;
                sw(o, &n.population);
                sw(o, &n.time);
            }
        };

        template <typename ostreamtype, std::size_t version>
        struct serialize_edge
        {
            inline void
            operator()(ostreamtype& o, const edge& e) const
            {
                throw std::runtime_error("invalid serialization version");
            }
        };

        template <typename ostreamtype>
        struct serialize_edge<ostreamtype, TS_TABLES_VERSION>
        {
            inline void
            operator()(ostreamtype& o, const edge& e) const
            {
                io::scalar_writer sw;
                sw(o, &e.left);
                sw(o, &e.right);
                sw(o, &e.parent);
                sw(o, &e.child);
            }
        };
    } // namespace ts
} // namespace fwdpp

#endif
