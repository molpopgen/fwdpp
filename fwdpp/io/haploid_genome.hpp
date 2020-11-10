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
        template <typename T> struct serialize_haploid_genome
        /// \brief Serialize a haploid_genome
        ///
        /// Serialize a haploid_genome. Specialize this for
        /// your haploid_genome type, else a static_assert
        /// will fail.  See implementation of
        /// serialize_haploid_genome<haploid_genome> for example.
        {
            template <typename streamtype>
            inline void
            operator()(streamtype&, const T&) const
            {
                static_assert(meta::always_false<T>::value,
                              "serialize_haploid_genome not implemented for type");
            }
        };

        template <typename T> struct deserialize_haploid_genome
        /// \brief Deserialize a haploid_genome
        ///
        /// Deserialize a haploid_genome. Specialize this for
        /// your haploid_genome type, else a static_assert will
        /// fail.  See implementation of
        /// deserialize_haploid_genome<haploid_genome> for example.
        {
            template <typename streamtype>
            inline T
            operator()(streamtype&) const
            {
                static_assert(meta::always_false<T>::value,
                              "deserialize_haploid_genome not implemented for type");
            }
        };

        template <typename GenomeContainerType, typename ostreamtype>
        void
        write_haploid_genomes(ostreamtype& buffer,
                              const GenomeContainerType& haploid_genomes)
        /// \brief Serialize a container of haploid_genomes.
        ///
        /// Works via specialization of serialize_haploid_genome.
        {
            static_assert(traits::is_haploid_genome<
                              typename GenomeContainerType::value_type>::value,
                          "GenomeContainerType must be a container of haploid_genomes");
            std::size_t nhaploid_genomes = haploid_genomes.size();
            scalar_writer writer;
            writer(buffer, &nhaploid_genomes);
            serialize_haploid_genome<typename GenomeContainerType::value_type>
                haploid_genome_writer;
            for (const auto& g : haploid_genomes)
                {
                    haploid_genome_writer(buffer, g);
                }
        }

        template <typename GenomeContainerType, typename istreamtype>
        void
        read_haploid_genomes(istreamtype& buffer, GenomeContainerType& haploid_genomes)
        /// \brief Deserialize a container of haploid_genomes.
        ///
        /// Works via specialization of deserialize_haploid_genome.
        {
            static_assert(traits::is_haploid_genome<
                              typename GenomeContainerType::value_type>::value,
                          "GenomeContainerType must be a container of haploid_genomes");
            std::size_t nhaploid_genomes;
            scalar_reader reader;
            reader(buffer, &nhaploid_genomes);
            deserialize_haploid_genome<typename GenomeContainerType::value_type>
                haploid_genome_reader;
            for (std::size_t i = 0; i < nhaploid_genomes; ++i)
                {
                    haploid_genomes.emplace_back(haploid_genome_reader(buffer));
                }
        }
    } // namespace io
} // namespace fwdpp
#endif
