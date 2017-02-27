#ifndef FWDPP_INTERNAL_TYPE_TRAITS_HPP
#define FWDPP_INTERNAL_TYPE_TRAITS_HPP
#include <type_traits>
#include <utility>

namespace KTfwd
{
    namespace traits
    {
        namespace internal
        {
            // Based on
            // http://stackoverflow.com/questions/11813940/possible-to-use-type-traits-sfinae-to-find-if-a-class-defines-a-member-typec

            template <typename...> struct void_t
            {
                typedef void type;
            };

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
/*
   This macro is used to generate
   types able to determine if a class type
   has a particular member type name
*/
#define HAS_TYPE(NAME)                                                        \
    template <typename, typename = void> struct has_##NAME : std::false_type  \
    {                                                                         \
    };                                                                        \
    template <typename T>                                                     \
    struct has_##NAME<T, typename void_t<typename T::NAME>::type>             \
        : std::true_type                                                      \
    {                                                                         \
    };

            HAS_TYPE(gamete_tag);
            HAS_TYPE(mutation_type);
            HAS_TYPE(first_type);
            HAS_TYPE(second_type);
        }
    }
}

#endif
