#ifndef FWDPP_IO_GAMETE_HPP__
#define FWDPP_IO_GAMETE_HPP__

#include <cstddef>
#include "scalar_serialization.hpp"

namespace fwdpp
{
    namespace io
    {
        template <typename T> struct serialize_gamete
        /// \brief Serialize a gamete
        ///
        /// Serialize a gamete. The implementation
		/// assumes fwdpp::gamete.  If you have derived
		/// a gamete from this type, then you must specialize
		/// this struct.
        {
            scalar_writer writer;
            serialize_gamete() : writer{} {}
            template <typename streamtype>
            inline void
            operator()(const T& g, streamtype& buffer) const
            {
                writer(buffer, &g.n);
                std::size_t nm = g.mutations.size();
                writer(buffer, &nm);
                if (nm)
                    {
                        writer(buffer, g.mutations.data(), nm);
                    }
                nm = g.smutations.size();
                writer(buffer, &nm);
                if (nm)
                    {
                        writer(buffer, g.smutations.data(), nm);
                    }
            }
        };

        template <typename T> struct deserialize_gamete
        /// \brief Deserialize a gamete
        ///
        /// Deserialize a gamete. The implementation
		/// assumes fwdpp::gamete.  If you have derived
		/// a gamete from this type, then you must specialize
		/// this struct.
        {
            scalar_reader reader;
            deserialize_gamete() : reader{} {}
            template <typename streamtype>
            inline T
            operator()(streamtype& buffer) const
            {
                decltype(T::n) n;
                std::size_t nm;
                decltype(T::mutations) mutations, smutations;
                reader(buffer, &n);
                reader(buffer, &nm);
                if (nm)
                    {
                        mutations.resize(nm);
                        reader(buffer, mutations.data(), nm);
                    }
                reader(buffer, &nm);
                if (nm)
                    {
                        smutations.resize(nm);
                        reader(buffer, smutations.data(), nm);
                    }
                return T(n, std::move(mutations), std::move(smutations));
            }
        };

        template <typename gcont_t, typename ostreamtype>
        void
        write_gametes(const gcont_t& gametes, ostreamtype& buffer)
        /// \brief Serialize a container of gametes.
        ///
        /// Works via argument-dependent lookup of serialize_gamete.
        {
            std::size_t ngametes = gametes.size();
            scalar_writer writer;
            writer(buffer, &ngametes);
            serialize_gamete<typename gcont_t::value_type> gamete_writer;
            for (const auto& g : gametes)
                {
                    gamete_writer(g, buffer);
                }
        }

        template <typename gcont_t, typename istreamtype>
        void
        read_gametes(gcont_t& gametes, istreamtype& buffer)
        /// \brief Deserialize a container of gametes.
        ///
        /// Works via argument-dependent lookup of deserialize_gamete.
        {
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
