#ifndef __FWDPP_INTERNAL_MUTATION_HPP__
#define __FWDPP_INTERNAL_MUTATION_HPP__

#include <algorithm>
#include <cassert>
namespace fwdpp
{
    namespace fwdpp_internal
    {
        template <typename mmodel, typename diploid_t, typename gamete_type,
                  typename mcont_type, typename queue_t>
        inline typename std::result_of<mmodel(queue_t &, mcont_type &)>::type
        mmodel_dispatcher(const mmodel &m, const diploid_t &,
                          const gamete_type &, mcont_type &mutations,
                          queue_t &recycling_bin)
        /*!
          Run-time dispatcher for mutation model
        */
        {
            return m(recycling_bin, mutations);
        }

        template <typename mmodel, typename diploid_t, typename gamete_type,
                  typename mcont_type, typename queue_t>
        inline typename std::result_of<mmodel(queue_t &, const gamete_type &,
                                              mcont_type &)>::type
        mmodel_dispatcher(const mmodel &m, const diploid_t &,
                          const gamete_type &g, mcont_type &mutations,
                          queue_t &recycling_bin)
        /*!
          Run-time dispatcher for mutation model
        */
        {
            return m(recycling_bin, g, mutations);
        }
        template <typename mmodel, typename diploid_t, typename gamete_type,
                  typename mcont_type, typename queue_t>
        inline typename std::result_of<mmodel(queue_t &, const diploid_t &,
                                              const gamete_type &,
                                              mcont_type &)>::type
        mmodel_dispatcher(const mmodel &m, const diploid_t &dip,
                          const gamete_type &g, mcont_type &mutations,
                          queue_t &recycling_bin)
        /*!
          Run-time dispatcher for mutation model
        */
        {
            return m(recycling_bin, dip, g, mutations);
        }
    }
}

#endif
