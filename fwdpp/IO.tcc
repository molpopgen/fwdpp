//  -*- C++ -*-
#ifndef __FWDPP_IO_TCC__
#define __FWDPP_IO_TCC__

#include <vector>
#include <algorithm>
#include <cassert>
#include <type_traits>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/io/scalar_serialization.hpp>
#include <fwdpp/io/diploid.hpp>
#include <fwdpp/internal/IOhelp.hpp>

namespace fwdpp
{
    // Binary I/O for individual-based simulation

    // Single-locus sims, single pop
    template <typename gcont_t, typename mcont_t, typename dipvector_t,
              typename ostreamtype>
    void
    write_binary_pop(const gcont_t &gametes, const mcont_t &mutations,
                     const dipvector_t &diploids, ostreamtype &buffer)
    {
        fwdpp_internal::write_mutations()(mutations, buffer);
        fwdpp_internal::write_haplotypes()(gametes, buffer);
        std::size_t NDIPS = diploids.size();
        io::scalar_writer writer;
        writer(buffer, &NDIPS);
		io::serialize_diploid<typename dipvector_t::value_type> dipwriter;
        for (const auto &dip : diploids)
            {
                dipwriter(dip, buffer);
            }
    }

    template <typename gcont_t, typename mcont_t, typename dipvector_t,
              typename istreamtype>
    void
    read_binary_pop(gcont_t &gametes, mcont_t &mutations,
                    dipvector_t &diploids, istreamtype &in)
    {
        gametes.clear();
        mutations.clear();
        diploids.clear();
        fwdpp_internal::read_mutations()(mutations, in);
        fwdpp_internal::read_haplotypes()(gametes, in);
        std::size_t NDIPS;
        io::scalar_reader()(in, &NDIPS);
        diploids.resize(NDIPS);
		io::deserialize_diploid<typename dipvector_t::value_type> dipreader;
        for (auto &dip : diploids)
            {
                dipreader(dip, in);
            }
    }

    // multi-locus, single pop, ostream
    template <typename gcont_t, typename mcont_t, typename dipvector_t,
              typename ostreamtype>
    void
    write_binary_pop_mloc(const gcont_t &mlocus_gametes,
                          const mcont_t &mutations,
                          const dipvector_t &diploids, ostreamtype &buffer)
    {
        unsigned nloci = unsigned(diploids[0].size());
        io::scalar_writer writer;
        writer(buffer, &nloci);
        // write mutations
        fwdpp_internal::write_mutations()(mutations, buffer);
        fwdpp_internal::write_haplotypes()(mlocus_gametes, buffer);
        unsigned ndips = unsigned(diploids.size());
        writer(buffer, &ndips);
		io::serialize_diploid<typename dipvector_t::value_type::value_type>
            dipwriter;
        for (const auto &dip : diploids)
            {
                for (const auto &genotype : dip)
                    {
                        dipwriter(genotype, buffer);
                    }
            }
    }

    // Multilocus, single-population, istream
    template <typename gcont_t, typename mcont_t, typename dipvector_t,
              typename istreamtype>
    void
    read_binary_pop_mloc(gcont_t &mlocus_gametes, mcont_t &mutations,
                         dipvector_t &diploids, istreamtype &in)
    {
        mlocus_gametes.clear();
        mutations.clear();
        diploids.clear();

        unsigned nloci;
        io::scalar_reader()(in, &nloci);
        // Read the mutations from the buffer
        fwdpp_internal::read_mutations()(mutations, in);
        fwdpp_internal::read_haplotypes()(mlocus_gametes, in);
        unsigned ndips;
        io::scalar_reader()(in, &ndips);
        diploids.resize(ndips, typename dipvector_t::value_type(nloci));
		io::deserialize_diploid<typename dipvector_t::value_type::value_type>
            dipreader;
        for (auto &dip : diploids)
            {
                assert(dip.size() == nloci);
                for (auto &genotype : dip)
                    {
                        dipreader(genotype, in);
                    }
            }
    }

    template <typename gcont_t, typename mcont_t, typename dipvector_t,
              typename ostreamtype>
    void
    write_binary_metapop(const gcont_t &gametes, const mcont_t &mutations,
                         const dipvector_t &diploids, ostreamtype &buffer)
    {
        std::size_t i = unsigned(diploids.size());
        io::scalar_writer writer;
        writer(buffer, &i);
        fwdpp_internal::write_mutations()(mutations, buffer);
        fwdpp_internal::write_haplotypes()(gametes, buffer);
		io::serialize_diploid<typename dipvector_t::value_type::value_type>
            dipwriter;
        for (const auto &deme : diploids)
            {
                i = deme.size();
                writer(buffer, &i);
                for (const auto &dip : deme)
                    {
                        dipwriter(dip, buffer);
                    }
            }
    }

    template <typename gcont_t, typename mcont_t, typename dipvector_t,
              typename istreamtype>
    void
    read_binary_metapop(gcont_t &gametes, mcont_t &mutations,
                        dipvector_t &diploids, istreamtype &in)
    {
        gametes.clear();
        mutations.clear();
        diploids.clear();

        std::size_t i;
        io::scalar_reader()(in, &i);
        diploids.resize(i);
        fwdpp_internal::read_mutations()(mutations, in);
        fwdpp_internal::read_haplotypes()(gametes, in);
		io::deserialize_diploid<typename dipvector_t::value_type::value_type>
            dipreader;
        for (auto &deme : diploids)
            {
                io::scalar_reader()(in, &i);
                if (i)
                    {
                        deme.resize(i);
                        for (auto &dip : deme)
                            {
                                dipreader(dip, in);
                            }
                    }
            }
    }

} // ns fwdpp

#endif
