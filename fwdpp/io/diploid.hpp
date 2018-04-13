#ifndef FWDPP_IO_DIPLOID_HPP__
#define FWDPP_IO_DIPLOID_HPP__

#include <cstdint>
#include <fwdpp/type_traits.hpp>
#include <fwdpp/meta/always_false.hpp>
#include "scalar_serialization.hpp"

// This file relies on std::uint32_t == fwdpp::uint_t,
// which prevents us from having to include that header.

namespace fwdpp
{
    namespace io
    {
        template <typename T> struct serialize_diploid
        /// \brief Serialize a diploid
        ///
        /// Serialize a diploid.
        /// This type simply writes T::first and T::second
        /// to a stream.  If your diploid type contains more data,
        /// then specialize this struct.
        /// \example custom_diploid.cc
        {
            template <typename ostreamtype>
            inline void
            operator()(ostreamtype &, const T &) const
            {
                static_assert(meta::always_false<T>::value,
                              "fwdpp::io::serialize_diploid not implemented "
                              "for this type");
            }
        };

        template <>
        struct serialize_diploid<std::pair<std::uint32_t, std::uint32_t>>
        /// \brief Serialize a diploid
        ///
        /// Serialize a diploid type std::pair<std::uint32_t, std::uint32_t>
        {
            io::scalar_writer writer;
            serialize_diploid() : writer{} {}
            template <typename ostreamtype>
            inline void
            operator()(
                ostreamtype &buffer,
                const std::pair<std::uint32_t, std::uint32_t> &dip) const
            {
                writer(buffer, &dip.first);
                writer(buffer, &dip.second);
            }
        };

        template <>
        struct serialize_diploid<std::pair<std::size_t, std::size_t>>
        /// \brief Serialize a diploid
        ///
        /// Serialize a diploid type std::pair<std::size_t, std::size_t>
        {
            io::scalar_writer writer;
            serialize_diploid() : writer{} {}
            template <typename ostreamtype>
            inline void
            operator()(ostreamtype &buffer,
                       const std::pair<std::size_t, std::size_t> &dip) const
            {
                writer(buffer, &dip.first);
                writer(buffer, &dip.second);
            }
        };

        template <typename T> struct deserialize_diploid
        /// \brief Deserialize a diploid
        ///
        /// Deserialize a diploid.
        /// \example custom_diploid.cc
        {
            template <typename istreamtype>
            inline void
            operator()(istreamtype &, T &) const
            {
                static_assert(meta::always_false<T>::value,
                              "fwdpp::io::deserialize_diploid not implemented "
                              "for this type");
            }
        };

        template <>
        struct deserialize_diploid<std::pair<std::uint32_t, std::uint32_t>>
        /// \brief Deserialize a diploid
        ///
        /// Deserialize a diploid type std::pair<std::uint32_t,std::uint32_t>
        {
            template <typename istreamtype>
            inline void
            operator()(istreamtype &buffer,
                       std::pair<std::uint32_t, std::uint32_t> &dip) const
            {
                io::scalar_reader reader;
                std::uint32_t c;
                reader(buffer, &c);
                dip.first = c;
                reader(buffer, &c);
                dip.second = c;
            }
        };

        template <>
        struct deserialize_diploid<std::pair<std::size_t, std::size_t>>
        /// \brief Deserialize a diploid
        ///
        /// Deserialize a diploid type std::pair<std::size_t,std::size_t>
        {
            template <typename istreamtype>
            inline void
            operator()(istreamtype &buffer,
                       std::pair<std::size_t, std::size_t> &dip) const
            {
                io::scalar_reader reader;
                std::size_t c;
                reader(buffer, &c);
                dip.first = c;
                reader(buffer, &c);
                dip.second = c;
            }
        };

        template <typename dipvector_t, typename ostreamtype>
        void
        write_diploids(ostreamtype &buffer, const dipvector_t &diploids)
        /// \brief Serialize a container of diploids.
        ///
        /// Works via specialization of serialize_diploid.
        {
            static_assert(
                traits::is_diploid<typename dipvector_t::value_type>::value,
                "dipvector_t must be a container of diploids");
            std::size_t NDIPS = diploids.size();
            io::scalar_writer writer;
            writer(buffer, &NDIPS);
            io::serialize_diploid<typename dipvector_t::value_type> dipwriter;
            for (const auto &dip : diploids)
                {
                    dipwriter(buffer, dip);
                }
        }

        template <typename dipvector_t, typename istreamtype>
        void
        read_diploids(istreamtype &buffer, dipvector_t &diploids)
        /// \brief Deserialize a container of diploids.
        ///
        /// Works via specialization of deserialize_diploid.
        {
            static_assert(
                traits::is_diploid<typename dipvector_t::value_type>::value,
                "dipvector_t must be a container of diploids");
            std::size_t NDIPS;
            io::scalar_reader reader;
            reader(buffer, &NDIPS);
            diploids.resize(NDIPS);
            io::deserialize_diploid<typename dipvector_t::value_type>
                dipreader;
            for (auto &dip : diploids)
                {
                    dipreader(buffer, dip);
                }
        }
    }
}

#endif
