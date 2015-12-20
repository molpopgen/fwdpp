#ifndef __FWDPP_TYPE_TRAITS_HPP__
#define __FWDPP_TYPE_TRAITS_HPP__

#include <fwdpp/internal/type_traits.hpp>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/tags/diploid_tags.hpp>
namespace KTfwd {
  namespace traits {
    //! Evaluates to std::true_type if T publicly inherits from KTfwd::tags::custom_diploid_t
    template<typename T>
    using is_custom_diploid_t = typename std::is_base_of<KTfwd::tags::custom_diploid_t,T>::type; 
  }
} 
#endif
