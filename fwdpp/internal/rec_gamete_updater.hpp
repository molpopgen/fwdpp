#ifndef __FWDPP_INTERNAL_REC_GAMETE_UPDATER_HPP__
#define __FWDPP_INTERNAL_REC_GAMETE_UPDATER_HPP__

#include <algorithm>
#include <functional>
#include <cassert>
namespace KTfwd
{
    namespace fwdpp_internal
    {
        template <typename itr_type, typename mcont_t>
        inline itr_type
        rec_update_itr(itr_type __first, itr_type __last,
                       const mcont_t &mutations, const double &val)
        {
            if (__first == __last)
                return __first;
            return std::upper_bound(
                __first, __last, std::cref(val),
                [&mutations](const double __val,
                             const std::size_t __mut) noexcept {
                    assert(__mut < mutations.size());
                    return __val < mutations[__mut].pos;
                });
        }

        template <typename itr_type, typename mcont_t,
                  typename mutation_index_cont_t>
        itr_type
        rec_gam_updater(itr_type __first, itr_type __last,
                        const mcont_t &mutations, mutation_index_cont_t &muts,
                        const double &val)
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
    }
}

#endif
