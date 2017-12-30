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

        struct write_haplotypes
        {
            template <typename gcont_t, typename ostreamtype>
            void
            operator()(const gcont_t &gametes, ostreamtype &buffer) const
            {
                std::size_t N = gametes.size();
                fwdpp::io::scalar_writer writer;
                writer(buffer, &N);
                for (const auto &g : gametes)
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

        struct read_haplotypes
        {
            template <typename gcont_t, typename istreamtype>
            void
            operator()(gcont_t &gametes, istreamtype &in) const
            {
                fwdpp::io::scalar_reader reader;
                std::size_t NHAPS;
                reader(in, &NHAPS);
                uint_t N;
                std::size_t nm;
                for (uint_t i = 0; i < NHAPS; ++i)
                    {
                        reader(in, &N);
                        typename gcont_t::value_type g(N);
                        reader(in, &nm);
                        if (nm)
                            {
                                g.mutations.resize(nm);
                                reader(in, g.mutations.data(), nm);
                            }
                        reader(in, &nm);
                        if (nm)
                            {
                                g.smutations.resize(nm);
                                reader(in, g.smutations.data(), nm);
                            }
                        gametes.emplace_back(std::move(g));
                    }
            }
        };
    }
}

#endif
