#ifndef __FWDPP_INTERNAL_TYPE_TRAITS_HPP__
#define __FWDPP_INTERNAL_TYPE_TRAITS_HPP__

#include <type_traits>

namespace KTfwd {
  namespace traits {
    namespace internal {

      /*
	Details of how policy arity can be checked at compile-time

	Based on
	http://stackoverflow.com/questions/20768649/standard-method-for-determining-the-arity-and-other-traits-of-stdbind-result,
	but rearranged to inherit from the C++11 integral constant, which make things a lot better when used in combination
	with std::enable_if, etc.
      */
      template<typename...> using void_t = void;

      template<typename F, typename enable = void>
      struct is_nullary : std::false_type {};
      
      template<typename F>
      struct is_nullary<F,void_t<typename std::result_of<F()>::type> > : std::true_type {};

      template<typename F, typename A1, typename enable = void>
      struct is_unary : std::false_type {};
      
      template<typename F, typename A1>
      struct is_unary<F,A1,void_t<typename std::result_of<F(A1)>::type> > : std::true_type {};

      template<typename F, typename A1, typename A2, typename enable = void>
      struct is_binary : std::false_type {};
      
      template<typename F, typename A1,typename A2>
      struct is_binary<F,A1,A2,void_t<typename std::result_of<F(A1,A2)>::type> > : std::true_type {};
      
    }
  }
}

#endif
