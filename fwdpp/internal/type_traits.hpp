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

            template <typename dipvector_t, typename gcont_t, typename mcont_t,
                      typename = void>
            struct fitness_fxn
            {
                using type = void;
            };

            template <typename dipvector_t, typename gcont_t, typename mcont_t>
            struct fitness_fxn<dipvector_t, gcont_t, mcont_t,
                               typename void_t<
                                   typename dipvector_t::value_type,
                                   typename gcont_t::value_type,
                                   typename mcont_t::value_type>::type>
            {
                using type = typename std::
                    conditional<(is_diploid<
                                     typename dipvector_t::value_type>::value
                                 || is_multilocus_diploid<
                                        typename dipvector_t::value_type>::
                                        value)
                                    && is_gamete<
                                           typename gcont_t::value_type>::value
                                    && is_mutation<typename mcont_t::
                                                       value_type>::value,
                                std::function<double(
                                    const typename dipvector_t::value_type &,
                                    const gcont_t &, const mcont_t &)>,
                                void>::type;
            };

        template <typename ff, typename dipvector_t, typename gcont_t,
                  typename mcont_t,
                  typename ffxn_t
                  = typename fitness_fxn<dipvector_t, gcont_t, mcont_t>::type>
        struct is_fitness_fxn
            : std::integral_constant<bool,
                                     !std::is_void<ffxn_t>::value
                                         && std::is_convertible<ff,
                                                                ffxn_t>::value>
        {
        };
        }
    }
}

#endif
