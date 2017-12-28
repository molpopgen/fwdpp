/*!
 * \file fwdpp/internal/type_traits.hpp
 * \not This file is not self-contained and cannot be included directly
 * Nasty SFINAE details of namespace KTfwd::traits
 */
#ifndef FWDPP_INTERNAL_TYPE_TRAITS_HPP
#define FWDPP_INTERNAL_TYPE_TRAITS_HPP
#include <type_traits>
#include <functional>
#include <utility>
#include <fwdpp/internal/void_t.hpp>

namespace KTfwd
{
    namespace traits
    {
        namespace internal
        {
            /* The SFINAE scheme will be documented by looking at
             * how traits::internal::is_diploid is implemented
             */

            template <typename T, typename = void>
            struct is_diploid : std::false_type
            /* Fallback type.  Evaluates to std::false_type */
            {
            };

            template <typename T>
            struct is_diploid<T, typename traits::internal::void_t<
                                     typename T::first_type,
                                     typename T::second_type>::type>
                /* If T does not have member type name first_type and
                 * second_type, this will fail to compile, and the fallback
                 * template type will be used.
                                 *
                                 * If successful, this type evaluates to
                 * std::true_type
                                 * if and only if first_type and second_type
                 * are the same
                                 * integer types, thus passing the minimal API
                 * requirement
                                 * for a diploid.
                 */
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

            template <
                typename ff, typename dipvector_t, typename gcont_t,
                typename mcont_t,
                typename ffxn_t =
                    typename fitness_fxn<dipvector_t, gcont_t, mcont_t>::type>
            struct is_fitness_fxn
                : std::integral_constant<bool,
                                         !std::is_void<ffxn_t>::value
                                             && std::is_convertible<ff,
                                                                    ffxn_t>::
                                                    value>
            {
            };

            template <typename mmodel_t, typename mcont_t, typename gcont_t,
                      typename = void, typename = void, typename = void>
            struct is_mutation_model : std::false_type
            {
            };

            template <typename mmodel_t, typename mcont_t, typename gcont_t>
            struct is_mutation_model<mmodel_t, mcont_t, gcont_t,
                                     typename void_t<
                                         typename std::result_of<mmodel_t(
                                             recycling_bin_t<mcont_t> &,
                                             mcont_t &)>::type>::type,
                                     typename std::enable_if<is_mutation<
                                         typename mcont_t::value_type>::
                                                                 value>::type,
                                     typename std::enable_if<is_gamete<
                                         typename gcont_t::value_type>::
                                                                 value>::type>
                : std::true_type
            {
            };
            template <typename mmodel_t, typename mcont_t, typename gcont_t>
            struct is_mutation_model<mmodel_t, mcont_t, gcont_t,
                                     typename void_t<
                                         typename std::result_of<mmodel_t(
                                             recycling_bin_t<mcont_t> &,
                                             typename gcont_t::value_type &,
                                             mcont_t &)>::type>::type,
                                     typename std::enable_if<is_mutation<
                                         typename mcont_t::value_type>::
                                                                 value>::type,
                                     typename std::enable_if<is_gamete<
                                         typename gcont_t::value_type>::
                                                                 value>::type>
                : std::true_type
            {
            };

            template <typename gcont_t_or_gamete_t, typename mcont_t,
                      typename = void>
            struct rich_recmodel_t
            {
                using type = typename std::
                    conditional<is_gamete<gcont_t_or_gamete_t>::value,
                                std::function<std::vector<double>(
                                    const gcont_t_or_gamete_t &,
                                    const gcont_t_or_gamete_t &,
                                    const mcont_t &)>,
                                void>::type;
            };

            template <typename gcont_t_or_gamete_t, typename mcont_t>
            struct rich_recmodel_t<gcont_t_or_gamete_t, mcont_t,
                                   typename void_t<
                                       typename gcont_t_or_gamete_t::
                                           value_type>::type>
            {
                using type = typename std::
                    conditional<is_gamete<typename gcont_t_or_gamete_t::
                                              value_type>::value,
                                std::function<std::vector<double>(
                                    const typename gcont_t_or_gamete_t::
                                        value_type &,
                                    const typename gcont_t_or_gamete_t::
                                        value_type &,
                                    const mcont_t &)>,
                                void>::type;
            };

            template <typename recmodel_t, typename diploid_t,
                      typename gamete_t, typename mcont_t, typename = void,
                      typename = void, typename = void, typename = void>
            struct is_rec_model : std::false_type
            {
            };

            template <typename recmodel_t, typename diploid_t,
                      typename gamete_t, typename mcont_t>
            struct is_rec_model<recmodel_t, diploid_t, gamete_t, mcont_t,
                                typename void_t<
                                    typename std::result_of<recmodel_t()>::
                                        type>::type,
                                typename std::enable_if<is_diploid<diploid_t>::
                                                            value>::type,
                                typename std::enable_if<is_gamete<gamete_t>::
                                                            value>::type,
                                typename std::enable_if<is_mutation<
                                    typename mcont_t::value_type>::value>::
                                    type> : std::true_type
            {
            };

            template <typename recmodel_t, typename diploid_t,
                      typename gamete_t, typename mcont_t>
            struct is_rec_model<recmodel_t, diploid_t, gamete_t, mcont_t,
                                typename void_t<
                                    typename std::result_of<recmodel_t(
                                        const gamete_t &, const gamete_t &,
                                        const mcont_t &)>::type>::type,
                                typename std::enable_if<is_diploid<diploid_t>::
                                                            value>::type,
                                typename std::enable_if<is_gamete<gamete_t>::
                                                            value>::type,
                                typename std::enable_if<is_mutation<
                                    typename mcont_t::value_type>::value>::
                                    type> : std::true_type
            {
            };

            template <typename recmodel_t, typename diploid_t,
                      typename gamete_t, typename mcont_t>
            struct is_rec_model<recmodel_t, diploid_t, gamete_t, mcont_t,
                                typename void_t<
                                    typename std::result_of<recmodel_t(
                                        const diploid_t &, const gamete_t &,
                                        const gamete_t &,
                                        const mcont_t &)>::type>::type,
                                typename std::enable_if<is_diploid<diploid_t>::
                                                            value>::type,
                                typename std::enable_if<is_gamete<gamete_t>::
                                                            value>::type,
                                typename std::enable_if<is_mutation<
                                    typename mcont_t::value_type>::value>::
                                    type> : std::true_type
            {
            };

            template <typename mcont_t, typename = void> struct mmodel_t
            {
                using type = void;
            };

            template <typename mcont_t>
            struct mmodel_t<mcont_t, typename void_t<
                                         typename mcont_t::value_type>::type>
            {
                using type = typename std::
                    conditional<is_mutation<
                                    typename mcont_t::value_type>::value,
                                std::function<std::size_t(
                                    recycling_bin_t<mcont_t> &, mcont_t &)>,
                                void>::type;
            };

            template <typename mcont_t, typename gcont_t, typename = void>
            struct mmodel_gamete_t
            {
                using type = void;
            };

            template <typename mcont_t, typename gcont_t>
            struct mmodel_gamete_t<mcont_t, gcont_t,
                                   typename void_t<
                                       typename mcont_t::value_type,
                                       typename gcont_t::value_type>::type>
            {
                using type = typename std::
                    conditional<is_mutation<
                                    typename mcont_t::value_type>::value
                                    && is_gamete<typename gcont_t::
                                                     value_type>::value,
                                std::function<std::size_t(
                                    recycling_bin_t<mcont_t> &,
                                    typename gcont_t::value_type &,
                                    mcont_t &)>,
                                void>::type;
            };
        }
    }
}

#endif
