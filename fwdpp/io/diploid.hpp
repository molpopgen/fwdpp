#ifndef FWDPP_IO_DIPLOID_HPP__
#define FWDPP_IO_DIPLOID_HPP__

#include <fwdpp/type_traits.hpp>
#include "scalar_serialization.hpp"

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
        {
            io::scalar_writer writer;
            serialize_diploid() : writer{} {}
            template <typename ostreamtype>
            inline void
            operator()(ostreamtype &buffer, const T &dip) const
            {
                writer(buffer, &dip.first);
                writer(buffer, &dip.second);
            }
        };

        template <typename T> struct deserialize_diploid
        /// \brief Deserialize a diploid
        ///
        /// Deserialize a diploid.
        /// This type simply reads T::first and T::second
        /// from a stream.  If your diploid type contains more data,
        /// then specialize this struct.
        {
            io::scalar_reader reader;
            deserialize_diploid() : reader{} {}
            template <typename istreamtype>
            inline void
            operator()(istreamtype &buffer, T &dip) const
            {
                typename T::first_type c;
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
