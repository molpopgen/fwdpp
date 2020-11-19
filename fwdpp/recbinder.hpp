#ifndef FWDPP_RECBINDER_HPP__
#define FWDPP_RECBINDER_HPP__

#include <vector>
#include <functional>
#include <gsl/gsl_rng.h>

namespace fwdpp
{
    template <typename T>
    inline std::function<std::vector<double>(void)>
    recbinder(T&& recmodel, const gsl_rng* r)
    /*! Convenience utility for binding simple recombination models.
     *  \version 0.6.0
     *  First added to library
     */
    {
        return [recmodel, r]() { return recmodel(r); };
    }
}

#endif
