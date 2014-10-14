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
  /*! \brief Returns true if member n == 0
    Used internally
   */
  struct n_is_zero
  {
    typedef bool result_type;
    template<typename T>
    inline bool operator()(const T & t)const
    {
      return t.n==0;
    }
  };

  /*! \brief Returns true if std::fabs(m.pos-d) <= std::numeric_limits<double>::epsilon()
    Returns true if std::fabs(m.pos-d) <= std::numeric_limits<double>::epsilon()
   */
  struct mutation_at_pos
  {
    typedef bool result_type;
    template<typename mutation_type>
    inline bool operator()(const mutation_type & m, const double & d) const
    {
      return( std::fabs(m.pos-d) <= std::numeric_limits<double>::epsilon() );
    }
  };
  
  /*!\brief Returns true if std::fabs(m1.pos-m2.pos) <= std::numeric_limits<double>::epsilon()
    Returns true if std::fabs(m1.pos-m2.pos) <= std::numeric_limits<double>::epsilon()
  */
  struct same_pos
  {
    typedef bool result_type;
    template<typename mutation_type>
    inline bool operator()(const mutation_type & m1, const mutation_type & m2) const
    {
      return( std::fabs(m1.pos-m2.pos) <= std::numeric_limits<double>::epsilon() );
    }
  };

  /// \brief Returns true if std::max(lhs,rhs)-std::min(lhs,rhs) <= std::numeric_limits<T>::epsilon()
  struct equal_eps
  {
    /*! \brief Returns true if std::max(lhs,rhs)-std::min(lhs,rhs) <= std::numeric_limits<T>::epsilon()
      Returns true if std::max(lhs,rhs)-std::min(lhs,rhs) <= std::numeric_limits<T>::epsilon()
    */
    typedef bool result_type;
    template<typename T>
    inline bool operator()(const T & lhs, const T & rhs) const
    {
      return( std::max(lhs,rhs)-std::min(lhs,rhs) <= std::numeric_limits<T>::epsilon() );
    }
  };

  /*! \brief Policy determining how mutations are removed from gametes after sampling
    Policy determining how mutations are removed from gametes after sampling

    Typical usage is to pass to KTfwd::sample_diploid as follows:
    \code
    std::bind( KTfwd::mutation_remover(), 0 ); //If a mutation is at frequency zero after sampling, remove pointers to it from all gametes
    std::bind( KWfwd::mutation_remover(), 0 ,twoN); //If a mutation is at frequency zero or 2N after sampling, remove pointers to it from all gametes
    \endcode
   */
  struct mutation_remover
  {
    typedef bool result_type;
    template<typename iterator_type> 
    inline result_type operator()(const iterator_type & i,
				  const unsigned & x1 ) const
    {
      return i->n == x1;
    }
    template<typename iterator_type>
    inline result_type operator()(const iterator_type & i,
				  const unsigned & x1,
				  const unsigned & x2) const
    {
      return i->n == x1 || i->n == x2;
    }
  };
}
#endif /* _FWD_FUNCTIONAL_HPP_ */
