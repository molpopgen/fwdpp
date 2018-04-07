#ifndef FWDPP_FORWARD_TYPES_SERIALIZATION_HPP__
#define FWDPP_FORWARD_TYPES_SERIALIZATION_HPP__

#include <fwdpp/io/mutation.hpp>
#include <fwdpp/io/gamete.hpp>

namespace fwdpp
{
    namespace io
    {
        template <> struct serialize_mutation<mutation>
        {
            template <typename streamtype>
            inline void
            operator()(streamtype &buffer, const mutation &m) const
            {
                io::scalar_writer writer;
                writer(buffer, &m.pos);
                writer(buffer, &m.s);
                writer(buffer, &m.h);
            }
        };

        template <> struct deserialize_mutation<mutation>
        {
            template <typename streamtype>
            inline mutation
            operator()(streamtype &buffer) const
            {
                double pos, s, h;
                io::scalar_reader reader;
                reader(buffer, &pos);
                reader(buffer, &s);
                reader(buffer, &h);

                return mutation(pos, s, h);
            }
        };
        template <> struct serialize_gamete<gamete>
        /// \brief Serialize a gamete
        ///
        /// Serialize fwdpp::gamete
        {
            template <typename streamtype>
            inline void
            operator()(streamtype &buffer, const gamete &g) const
            {
                scalar_writer writer;
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
        template <> struct deserialize_gamete<gamete>
        /// \brief Deserialize a gamete
        ///
        /// Deserialize a fwdpp::gamete.
        {
            template <typename streamtype>
            inline gamete
            operator()(streamtype& buffer) const
            {
                scalar_reader reader;
                decltype(gamete::n) n;
                std::size_t nm;
                decltype(gamete::mutations) mutations, smutations;
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
                return gamete(n, std::move(mutations), std::move(smutations));
            }
        };
    }
}

#endif
