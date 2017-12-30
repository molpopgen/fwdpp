#ifndef FWDPP_IO_MUTATION_HPP__
#define FWDPP_IO_MUTATION_HPP__

#include <cstddef>
#include <stdexcept>
#include <fwdpp/forward_types.hpp>
#include "scalar_serialization.hpp"

/*! \namespace fwdpp::io
 * \brief I/O types and functions
 */

namespace fwdpp
{
    namespace io
    {
        template <typename T> struct serialize_mutation
        {
            template <typename istreamtype>
            inline void
            operator()(const T &, istreamtype &) const
            {
                throw std::runtime_error(
                    "serializtion not implemented for this mutation type");
            }
        };

        template <typename T> struct deserialize_mutation
        {
            template <typename ostreamtype>
            inline T
            operator()(ostreamtype &) const
            {
                throw std::runtime_error(
                    "deserializtion not implemented for this mutation type");
            }
        };

        template <typename mcont_t, typename ostreamtype>
        void
        write_mutations(const mcont_t &mutations, ostreamtype &buffer)
        {
            std::size_t MUTNO = mutations.size();
            fwdpp::io::scalar_writer()(buffer, &MUTNO);
            // write the mutation data to the buffer
            fwdpp::io::serialize_mutation<typename mcont_t::value_type> mw;
            for (const auto &m : mutations)
                mw(m, buffer);
        }

        template <typename mcont_t, typename istreamtype>
        void
        read_mutations(mcont_t &mutations, istreamtype &in)
        {
            std::size_t NMUTS;
            fwdpp::io::scalar_reader()(in, &NMUTS);
            fwdpp::io::deserialize_mutation<typename mcont_t::value_type> mr;
            for (uint_t i = 0; i < NMUTS; ++i)
                {
                    mutations.emplace_back(mr(in));
                }
        }
    }
}

#endif
