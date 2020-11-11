#ifndef __FDWPP_INTERNAL_RECOMBINATION_COMMON_HPP__
#define __FDWPP_INTERNAL_RECOMBINATION_COMMON_HPP__

#include <fwdpp/debug.hpp>
#include <fwdpp/internal/rec_gamete_updater.hpp>
#include <cassert>
namespace fwdpp
{

    namespace fwdpp_internal
    {

        template <typename double_vec_type, typename GenomeContainerType,
                  typename MutationContainerType>
        void
        recombine_haploid_genomes(
            const double_vec_type &pos, const std::size_t ibeg, const std::size_t jbeg,
            const GenomeContainerType &haploid_genomes,
            const MutationContainerType &mutations,
            typename GenomeContainerType::value_type::mutation_container &neutral,
            typename GenomeContainerType::value_type::mutation_container &selected)
        {
            assert(std::is_sorted(pos.cbegin(), pos.cend()));

            auto itr = haploid_genomes[ibeg].mutations.cbegin();
            auto jtr = haploid_genomes[jbeg].mutations.cbegin();
            auto itr_s = haploid_genomes[ibeg].smutations.cbegin();
            auto jtr_s = haploid_genomes[jbeg].smutations.cbegin();
            auto itr_e = haploid_genomes[ibeg].mutations.cend();
            auto itr_s_e = haploid_genomes[ibeg].smutations.cend();
            auto jtr_e = haploid_genomes[jbeg].mutations.cend();
            auto jtr_s_e = haploid_genomes[jbeg].smutations.cend();

            for (auto &&p : pos)
                {
                    itr = fwdpp_internal::rec_gam_updater(itr, itr_e, mutations, neutral,
                                                          p);
                    itr_s = fwdpp_internal::rec_gam_updater(itr_s, itr_s_e, mutations,
                                                            selected, p);
                    jtr = fwdpp_internal::rec_update_itr(jtr, jtr_e, mutations, p);
                    jtr_s = fwdpp_internal::rec_update_itr(jtr_s, jtr_s_e, mutations, p);
                    std::swap(itr, jtr);
                    std::swap(itr_s, jtr_s);
                    std::swap(itr_e, jtr_e);
                    std::swap(itr_s_e, jtr_s_e);
                }
            debug::haploid_genome_is_sorted(haploid_genomes[ibeg], mutations);
        }
    }
}
#endif
