#ifndef FWDPP_FORWARD_TYPES_SERIALIZATION_HPP__
#define FWDPP_FORWARD_TYPES_SERIALIZATION_HPP__

#include <fwdpp/io/mutation.hpp>
#include <fwdpp/io/haploid_genome.hpp>

namespace fwdpp
{
    namespace io
    {
        template <> struct serialize_haploid_genome<haploid_genome>
        /// \brief Serialize a haploid_genome
        ///
        /// Serialize fwdpp::haploid_genome
        {
            template <typename streamtype>
            inline void
            operator()(streamtype &buffer, const haploid_genome &g) const
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
        template <> struct deserialize_haploid_genome<haploid_genome>
        /// \brief Deserialize a haploid_genome
        ///
        /// Deserialize a fwdpp::haploid_genome.
        {
            template <typename streamtype>
            inline haploid_genome
            operator()(streamtype &buffer) const
            {
                scalar_reader reader;
                decltype(haploid_genome::n) n;
                std::size_t nm;
                decltype(haploid_genome::mutations) mutations, smutations;
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
                return haploid_genome(n, std::move(mutations), std::move(smutations));
            }
        };
    }
}

#endif
