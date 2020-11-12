#ifndef __FWDPP_INTERNAL_MUTATION_HPP__
#define __FWDPP_INTERNAL_MUTATION_HPP__

#include <algorithm>

namespace fwdpp
{
    namespace fwdpp_internal
    {
        template <typename mmodel, typename DiploidType, typename haploid_genome_type,
                  typename MutationContainerType, typename queue_t>
        inline typename std::result_of<mmodel(queue_t &, MutationContainerType &)>::type
        mmodel_dispatcher(const mmodel &m, const DiploidType &,
                          const haploid_genome_type &, MutationContainerType &mutations,
                          queue_t &recycling_bin)
        /*!
          Run-time dispatcher for mutation model
        */
        {
            return m(recycling_bin, mutations);
        }

        template <typename mmodel, typename DiploidType, typename haploid_genome_type,
                  typename MutationContainerType, typename queue_t>
        inline typename std::result_of<mmodel(const haploid_genome_type &, queue_t &,
                                              MutationContainerType &)>::type
        mmodel_dispatcher(const mmodel &m, const DiploidType &,
                          const haploid_genome_type &g, MutationContainerType &mutations,
                          queue_t &recycling_bin)
        /*!
          Run-time dispatcher for mutation model
        */
        {
            return m(g, recycling_bin, mutations);
        }

        template <typename mmodel, typename DiploidType, typename haploid_genome_type,
                  typename MutationContainerType, typename queue_t>
        inline typename std::result_of<mmodel(queue_t &, const DiploidType &,
                                              const haploid_genome_type &,
                                              MutationContainerType &)>::type
        mmodel_dispatcher(const mmodel &m, const DiploidType &dip,
                          const haploid_genome_type &g, MutationContainerType &mutations,
                          queue_t &recycling_bin)
        /*!
          Run-time dispatcher for mutation model
        */
        {
            return m(recycling_bin, dip, g, mutations);
        }
    } // namespace fwdpp_internal
} // namespace fwdpp

#endif
