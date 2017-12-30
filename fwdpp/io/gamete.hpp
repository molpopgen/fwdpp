#ifndef FWDPP_IO_GAMETE_HPP__
#define FWDPP_IO_GAMETE_HPP__

#include <cstddef>
#include "scalar_serialization.hpp"

namespace fwdpp
{
    namespace io
    {
        template <typename T> struct serlialize_gamete
        {
            scalar_writer writer;
            serlialize_gamete() : writer{} {}
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

        template <typename T> struct deserlialize_gamete
        {
            scalar_reader reader;
            deserlialize_gamete() : reader{} {}
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
        {
            std::size_t ngametes = gametes.size();
            scalar_writer writer;
            writer(buffer, &ngametes);
            serlialize_gamete<typename gcont_t::value_type> gamete_writer;
            for (const auto& g : gametes)
                {
                    gamete_writer(g, buffer);
                }
        }

        template <typename gcont_t, typename istreamtype>
        void
        read_gametes(gcont_t& gametes, istreamtype& buffer)
        {
            std::size_t ngametes;
            scalar_reader reader;
            reader(buffer, &ngametes);
            deserlialize_gamete<typename gcont_t::value_type> gamete_reader;
            for (std::size_t i = 0; i < ngametes; ++i)
                {
					gametes.emplace_back(gamete_reader(buffer));
                }
        }
    }
}
#endif
