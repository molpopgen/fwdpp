#ifndef FWDPP_INTERNAL_REC_GAMETE_UPDATER_HPP
#define FWDPP_INTERNAL_REC_GAMETE_UPDATER_HPP

#include <algorithm>
#include <functional>

namespace fwdpp
{
    namespace fwdpp_internal
    {
        template <typename itr_type, typename mcont_t>
        inline itr_type
        rec_update_itr(itr_type __first, itr_type __last,
                       const mcont_t &mutations, const double val)
        {
            if (__first == __last)
                return __first;
            return std::lower_bound(
                __first, __last, val,
                [&mutations](const typename itr_type::value_type __mut,
                             const double v) {
                    return mutations[__mut].pos < v;
                });
        }

        template <typename itr_type, typename mcont_t,
                  typename mutation_index_cont_t>
        itr_type
        rec_gam_updater(itr_type __first, itr_type __last,
                        const mcont_t &mutations, mutation_index_cont_t &muts,
                        const double val)
        {
            // O(log_2) comparisons of double plus at most __last - __first
            // copies
            itr_type __ub = rec_update_itr(__first, __last, mutations, val);
            /*
              NOTE: the use of insert here
              instead of std::copy(__first,__ub,std::back_inserter(muts));
              Reduced peak RAM use on GCC 4.9.2/Ubuntu Linux.
            */
            muts.insert(muts.end(), __first, __ub);
            return __ub;
        }
    } // namespace fwdpp_internal
} // namespace fwdpp

#endif
