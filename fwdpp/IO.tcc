//  -*- C++ -*-
#ifndef __FWDPP_IO_TCC__
#define __FWDPP_IO_TCC__

#include <vector>
#include <algorithm>
#include <cassert>
#include <type_traits>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/internal/IOhelp.hpp>

namespace KTfwd
{
    // Binary I/O for individual-based simulation

    // Single-locus sims, single pop
    template <typename gcont_t, typename mcont_t, typename dipvector_t,
              typename mutation_writer_type, typename ostreamtype,
              typename diploid_writer_t>
    void
    write_binary_pop(const gcont_t &gametes, const mcont_t &mutations,
                     const dipvector_t &diploids,
                     const mutation_writer_type &mw, ostreamtype &buffer,
                     const diploid_writer_t &dw)
    {
        fwdpp_internal::write_mutations()(mutations, mw, buffer);
        fwdpp_internal::write_haplotypes()(gametes, buffer);
        std::size_t NDIPS = diploids.size();
        fwdpp_internal::scalar_writer writer;
        writer(buffer, &NDIPS);
        for (const auto &dip : diploids)
            {
                writer(buffer, &dip.first);
                writer(buffer, &dip.second);
                dw(dip, buffer);
            }
    }

    template <typename gcont_t, typename mcont_t, typename dipvector_t,
              typename mutation_reader_type, typename istreamtype,
              typename diploid_reader_t>
    void
    read_binary_pop(gcont_t &gametes, mcont_t &mutations,
                    dipvector_t &diploids, const mutation_reader_type &mr,
                    istreamtype &in, const diploid_reader_t &dr)
    {
        gametes.clear();
        mutations.clear();
        diploids.clear();
        fwdpp_internal::read_mutations()(mutations, mr, in);
        fwdpp_internal::read_haplotypes()(gametes, in);
        std::size_t NDIPS, c;
        fwdpp_internal::scalar_reader()(in, &NDIPS);
        diploids.resize(NDIPS);
        for (auto &dip : diploids)
            {
                fwdpp_internal::scalar_reader()(in, &c);
                dip.first = c;
                fwdpp_internal::scalar_reader()(in, &c);
                dip.second = c;
                dr(dip, in);
            }
    }

    // multi-locus, single pop, ostream
    template <typename gcont_t, typename mcont_t, typename dipvector_t,
              typename mutation_writer_type, typename ostreamtype,
              typename diploid_writer_t>
    void
    write_binary_pop_mloc(const gcont_t &mlocus_gametes,
                          const mcont_t &mutations,
                          const dipvector_t &diploids,
                          const mutation_writer_type &mw, ostreamtype &buffer,
                          const diploid_writer_t &dw)
    {
        unsigned nloci = unsigned(diploids[0].size());
        fwdpp_internal::scalar_writer writer;
        writer(buffer, &nloci);
        // write mutations
        fwdpp_internal::write_mutations()(mutations, mw, buffer);
        fwdpp_internal::write_haplotypes()(mlocus_gametes, buffer);
        unsigned ndips = unsigned(diploids.size());
        writer(buffer, &ndips);
        for (const auto &dip : diploids)
            {
                for (const auto &genotype : dip)
                    {
                        writer(buffer, &genotype.first);
                        writer(buffer, &genotype.second);
                        dw(genotype, buffer);
                    }
            }
    }

    // Multilocus, single-population, istream
    template <typename gcont_t, typename mcont_t, typename dipvector_t,
              typename mutation_reader_type, typename istreamtype,
              typename diploid_reader_t>
    void
    read_binary_pop_mloc(gcont_t &mlocus_gametes, mcont_t &mutations,
                         dipvector_t &diploids, const mutation_reader_type &mr,
                         istreamtype &in, const diploid_reader_t &dr)
    {
        mlocus_gametes.clear();
        mutations.clear();
        diploids.clear();

        unsigned nloci;
        fwdpp_internal::scalar_reader()(in, &nloci);
        // Read the mutations from the buffer
        fwdpp_internal::read_mutations()(mutations, mr, in);
        fwdpp_internal::read_haplotypes()(mlocus_gametes, in);
        unsigned ndips;
        fwdpp_internal::scalar_reader()(in, &ndips);
        diploids.resize(ndips, typename dipvector_t::value_type(nloci));
        for (auto &dip : diploids)
            {
                assert(dip.size() == nloci);
                for (auto &genotype : dip)
                    {
                        fwdpp_internal::scalar_reader()(in, &genotype.first);
                        fwdpp_internal::scalar_reader()(in, &genotype.second);
                        dr(genotype, in);
                    }
            }
    }

    template <typename gcont_t, typename mcont_t, typename dipvector_t,
              typename mutation_writer_type, typename ostreamtype,
              typename diploid_writer_t>
    void
    write_binary_metapop(const gcont_t &gametes, const mcont_t &mutations,
                         const dipvector_t &diploids,
                         const mutation_writer_type &mw, ostreamtype &buffer,
                         const diploid_writer_t &dw)
    {
        std::size_t i = unsigned(diploids.size());
        fwdpp_internal::scalar_writer writer;
        writer(buffer, &i);
        fwdpp_internal::write_mutations()(mutations, mw, buffer);
        fwdpp_internal::write_haplotypes()(gametes, buffer);
        for (const auto &deme : diploids)
            {
                i = deme.size();
                writer(buffer, &i);
                for (const auto &dip : deme)
                    {
                        writer(buffer, &dip.first);
                        writer(buffer, &dip.second);
                        dw(dip, buffer);
                    }
            }
    }

    template <typename gcont_t, typename mcont_t, typename dipvector_t,
              typename mutation_reader_type, typename istreamtype,
              typename diploid_reader_t>
    void
    read_binary_metapop(gcont_t &gametes, mcont_t &mutations,
                        dipvector_t &diploids, const mutation_reader_type &mr,
                        istreamtype &in, const diploid_reader_t &dr)
    {
        gametes.clear();
        mutations.clear();
        diploids.clear();

        std::size_t i;
        fwdpp_internal::scalar_reader()(in, &i);
        diploids.resize(i);
        fwdpp_internal::read_mutations()(mutations, mr, in);
        fwdpp_internal::read_haplotypes()(gametes, in);
        for (auto &deme : diploids)
            {
                fwdpp_internal::scalar_reader()(in, &i);
                if (i)
                    {
                        deme.resize(i);
                        for (auto &dip : deme)
                            {
                                fwdpp_internal::scalar_reader()(in,
                                                                &dip.first);
                                fwdpp_internal::scalar_reader()(in,
                                                                &dip.second);
                                dr(dip, in);
                            }
                    }
            }
    }

} // ns KTfwd

#endif
