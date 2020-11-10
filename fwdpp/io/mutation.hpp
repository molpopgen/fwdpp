#ifndef FWDPP_IO_MUTATION_HPP__
#define FWDPP_IO_MUTATION_HPP__

#include <cstddef>
#include <typeinfo>
#include <fwdpp/type_traits.hpp>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/meta/always_false.hpp>
#include "scalar_serialization.hpp"

/*! \namespace fwdpp::io
 * \brief I/O types and functions.
 *
 * Template typenames referring to streams assume input
 * or output streams compatible with those from namespace
 * std.
 *
 * See @ref md_md_serialization for details.
 */

namespace fwdpp
{
    namespace io
    {

        template <typename T> struct serialize_mutation
        /// \brief Serialize a mutation
        ///
        /// Serialize a mutation.
        /// This *must* be specialized for a specific mutation type,
        /// otherwise a static_assertion will fail.
        /// The type T must be a valid mutation type.
        /// See serialize_mutation<popgenmut> for example implementation
        /// \example custom_mutation.cc
        {
            template <typename ostreamtype>
            inline void
            operator()(ostreamtype &, const T &) const
            {
                static_assert(meta::always_false<T>::value,
                              "fwdpp::io::serialize_mutation "
                              "not implemented for this "
                              "type");
            }
        };

        template <typename T> struct deserialize_mutation
        /// \brief Deserialize a mutation
        ///
        /// Deserialize a mutation.
        /// This *must* be specialized for a specific mutation type,
        /// otherwise a static_assertion will fail.
        /// The type T must be a valid mutation type.
        /// See deserialize_mutation<popgenmut> for example implementation
        /// \example custom_mutation.cc
        {
            template <typename istreamtype>
            inline T
            operator()(istreamtype &) const
            {
                static_assert(meta::always_false<T>::value,
                              "fwdpp::io::deserialize_mutation "
                              "not implemented for this "
                              "type");
            }
        };

        template <typename MutationContainerType, typename ostreamtype>
        void
        write_mutations(ostreamtype &buffer, const MutationContainerType &mutations)
        /// \brief Serialize a container of mutations.
        ///
        /// Works via specialization of serialize_mutation.
        {
            static_assert(
                traits::is_mutation_v<typename MutationContainerType::value_type>,
                "MutationContainerType must be a container of mutations");
            std::size_t MUTNO = mutations.size();
            fwdpp::io::scalar_writer()(buffer, &MUTNO);
            // write the mutation data to the buffer
            fwdpp::io::serialize_mutation<typename MutationContainerType::value_type> mw;
            for (const auto &m : mutations)
                mw(buffer, m);
        }

        template <typename MutationContainerType, typename istreamtype>
        void
        read_mutations(istreamtype &in, MutationContainerType &mutations)
        /// \brief Deserialize a container of mutations.
        ///
        /// Works via specialization of deserialize_mutation.
        {
            static_assert(
                traits::is_mutation_v<typename MutationContainerType::value_type>,
                "MutationContainerType must be a container of mutations");
            std::size_t NMUTS;
            fwdpp::io::scalar_reader()(in, &NMUTS);
            fwdpp::io::deserialize_mutation<typename MutationContainerType::value_type> mr;
            for (uint_t i = 0; i < NMUTS; ++i)
                {
                    mutations.emplace_back(mr(in));
                }
        }
    }
}

#endif
