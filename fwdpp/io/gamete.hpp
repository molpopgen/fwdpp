#ifndef FWDPP_IO_GAMETE_HPP__
#define FWDPP_IO_GAMETE_HPP__

#include <cstddef>
#include <fwdpp/type_traits.hpp>
#include <fwdpp/meta/always_false.hpp>
#include "scalar_serialization.hpp"

namespace fwdpp
{
    namespace io
    {
        template <typename T> struct serialize_gamete
        /// \brief Serialize a gamete
        ///
        /// Serialize a gamete. Specialize this for
        /// your gamete type, else a static_assert
        /// will fail.  See implementation of
        /// serialize_gamete<gamete> for example.
        {
            template <typename streamtype>
            inline void
            operator()(streamtype&, const T&) const
            {
                static_assert(meta::always_false<T>::value,
                              "serialize_gamete not implemented for type");
            }
        };

        template <typename T> struct deserialize_gamete
        /// \brief Deserialize a gamete
        ///
        /// Deserialize a gamete. Specialize this for
        /// your gamete type, else a static_assert will
        /// fail.  See implementation of
        /// deserialize_gamete<gamete> for example.
        {
            template <typename streamtype>
            inline T
            operator()(streamtype&) const
            {
                static_assert(meta::always_false<T>::value,
                              "deserialize_gamete not implemented for type");
            }
        };

        template <typename gcont_t, typename ostreamtype>
        void
        write_gametes(ostreamtype& buffer, const gcont_t& gametes)
        /// \brief Serialize a container of gametes.
        ///
        /// Works via specialization of serialize_gamete.
        {
            static_assert(
                traits::is_gamete<typename gcont_t::value_type>::value,
                "gcont_t must be a container of gametes");
            std::size_t ngametes = gametes.size();
            scalar_writer writer;
            writer(buffer, &ngametes);
            serialize_gamete<typename gcont_t::value_type> gamete_writer;
            for (const auto& g : gametes)
                {
                    gamete_writer(buffer, g);
                }
        }

        template <typename gcont_t, typename istreamtype>
        void
        read_gametes(istreamtype& buffer, gcont_t& gametes)
        /// \brief Deserialize a container of gametes.
        ///
        /// Works via specialization of deserialize_gamete.
        {
            static_assert(
                traits::is_gamete<typename gcont_t::value_type>::value,
                "gcont_t must be a container of gametes");
            std::size_t ngametes;
            scalar_reader reader;
            reader(buffer, &ngametes);
            deserialize_gamete<typename gcont_t::value_type> gamete_reader;
            for (std::size_t i = 0; i < ngametes; ++i)
                {
                    gametes.emplace_back(gamete_reader(buffer));
                }
        }
    }
}
#endif
