#ifndef _FWD_FUNCTIONAL_HPP_
#define _FWD_FUNCTIONAL_HPP_

#include <cmath>
#include <limits>
#include <functional>

/*! \file fwd_functional.hpp
  Defines several function objects used both internally and by library users
*/

namespace KTfwd
{
  /// \brief Returns true if std::max(lhs,rhs)-std::min(lhs,rhs) <= std::numeric_limits<T>::epsilon()
  struct equal_eps
  {
    /*! \brief Returns true if std::max(lhs,rhs)-std::min(lhs,rhs) <= std::numeric_limits<T>::epsilon()
      Returns true if std::max(lhs,rhs)-std::min(lhs,rhs) <= std::numeric_limits<T>::epsilon()
    */
    using result_type = bool;
    template<typename T>
    inline bool operator()(const T & lhs, const T & rhs) const
    {
      return( std::max(lhs,rhs)-std::min(lhs,rhs) <= std::numeric_limits<T>::epsilon() );
    }
  };
  
  /* \brief Policy telling library not to remove any mutations from gametes after sampling
     \note This is an empty struct that functions as a dispatch tag for library internals
   */
  struct remove_nothing
  {
  };
}
#endif /* _FWD_FUNCTIONAL_HPP_ */
