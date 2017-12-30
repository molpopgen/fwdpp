#ifndef __FWDPP_INTERNAL_IOHELP_HPP__
#define __FWDPP_INTERNAL_IOHELP_HPP__

/*
  Mechanics of data serialization
  The various write_binary_pop and read_binary pop
  functions rely on these implementations
*/
#include <zlib.h>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/io/mutation.hpp>
#include <fwdpp/io/scalar_serialization.hpp>

namespace fwdpp
{
    namespace fwdpp_internal
    {

        struct write_mutations
        {
            template <typename mcont_t, typename ostreamtype>
            void
            operator()(const mcont_t &mutations, ostreamtype &buffer) const
            {
                std::size_t MUTNO = mutations.size();
                fwdpp::io::scalar_writer()(buffer, &MUTNO);
                // write the mutation data to the buffer
                fwdpp::io::serialize_mutation<typename mcont_t::value_type> mw;
                for (const auto &m : mutations)
                    mw(m, buffer);
            }
        };

        struct read_mutations
        {
            template <typename mcont_t, typename istreamtype>
            void
            operator()(mcont_t &mutations, istreamtype &in) const
            {
                std::size_t NMUTS;
                fwdpp::io::scalar_reader()(in, &NMUTS);
                fwdpp::io::deserialize_mutation<typename mcont_t::value_type>
                    mr;
                for (uint_t i = 0; i < NMUTS; ++i)
                    {
                        mutations.emplace_back(mr(in));
                    }
            }
        };
    }
}

#endif
