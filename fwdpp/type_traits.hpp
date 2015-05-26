#ifndef __FWDPP_TYPE_TRAITS_HPP__
#define __FWDPP_TYPE_TRAITS_HPP__

#include <fwdpp/internal/type_traits.hpp>
#include <fwdpp/forward_types.hpp>
namespace KTfwd {
  namespace traits {

    //! Evaluates to std::true_type if F is a function taking no arguments
    template<typename F>
    using is_nullary_t = typename internal::is_nullary<F>::type;

    //! Evaluates to std::true_type if F is a function a single argument of type A1
    template<typename F,typename A1>
    using is_unary_t = typename internal::is_unary<F,A1>::type;
   
    /*! Evaluates to std::true_type if F is a function taking arguments of types A1 and A2
      \note Due to the vaguaries of std::bind, distinguishing a binary function from a unary one 
      requires something like: is_binary_t<F,A1,A2>::value && !(is_unary_t<F,A1>::value||is_unary_t<F,A2>::value)
     */
    template<typename F,typename A1,typename A2>
    using is_binary_t = typename internal::is_binary<F,A1,A2>::type;

  }
} 
#endif
