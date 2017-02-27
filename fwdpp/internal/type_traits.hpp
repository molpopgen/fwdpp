#ifndef FWDPP_INTERNAL_TYPE_TRAITS_HPP
#define FWDPP_INTERNAL_TYPE_TRAITS_HPP
#include <type_traits>
#include <fwdpp/internal/void_t.hpp>
#include <utility>

namespace KTfwd
{
    namespace traits
    {
        namespace internal
        {
            template <typename T, typename = void>
            struct is_diploid : std::false_type
            {
            };

            template <typename T>
            struct is_diploid<T, typename traits::internal::void_t<
                                     typename T::first_type,
                                     typename T::second_type>::type>
                : std::integral_constant<bool,
                                         std::is_integral<
                                             typename T::first_type>::value
                                             && std::is_same<
                                                    typename T::first_type,
                                                    typename T::second_type>::
                                                    value>
            {
            };

            template <typename T, typename = void>
            struct is_multilocus_diploid : std::false_type
            {
            };

            template <typename T>
            struct is_multilocus_diploid<T, typename void_t<
                                                typename T::value_type>::type>
                : std::integral_constant<bool,
                                         is_diploid<
                                             typename T::value_type>::value>
            {
            };

            template <typename T, typename = void>
            struct is_custom_diploid : std::false_type
            {
            };

            template <typename T>
            struct is_custom_diploid<T, typename void_t<
                                            typename T::first_type,
                                            typename T::second_type>::type>
                : std::
                      integral_constant<bool,
                                        is_diploid<T>::value
                                            && !std::
                                                   is_same<std::pair<
                                                               typename T::
                                                                   first_type,
                                                               typename T::
                                                                   second_type>,
                                                           T>::value>
            {
            };
        }
    }
}

#endif
