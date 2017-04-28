#ifndef _FWD_FUNCTIONAL_HPP_
#define _FWD_FUNCTIONAL_HPP_

#include <cmath>
#include <limits>
#include <functional>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/*! \file fwd_functional.hpp
  Defines several function objects used both internally and by library users
*/

namespace KTfwd
{
    /// \brief Returns true if std::max(lhs,rhs)-std::min(lhs,rhs) <=
    /// std::numeric_limits<T>::epsilon()
    struct equal_eps
    {
        /*! \brief Returns true if std::max(lhs,rhs)-std::min(lhs,rhs) <=
          std::numeric_limits<T>::epsilon()
          Returns true if std::max(lhs,rhs)-std::min(lhs,rhs) <=
          std::numeric_limits<T>::epsilon()
        */
        using result_type = bool;
        template <typename T>
        inline bool
        operator()(const T &lhs, const T &rhs) const
        {
            return (std::max(lhs, rhs) - std::min(lhs, rhs)
                    <= std::numeric_limits<T>::epsilon());
        }
    };

    /* \brief Policy telling library not to remove any mutations from gametes
       after sampling
       \note This is an empty struct that functions as a dispatch tag for
       library internals
     */
    struct remove_nothing
    {
    };

    /*! \brief Policy telling library to remove neutral mutations from gametes
     * after sampling
     */
    struct remove_neutral
    {
        template <typename mtype>
        inline bool
        operator()(const mtype &m) const
        {
            return m.neutral;
        }
    };
}
#endif /* _FWD_FUNCTIONAL_HPP_ */
