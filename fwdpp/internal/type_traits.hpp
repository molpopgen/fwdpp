/*!
 * \file fwdpp/internal/type_traits.hpp
 * \not This file is not self-contained and cannot be included directly
 * Nasty SFINAE details of namespace fwdpp::traits
 */
#ifndef FWDPP_INTERNAL_TYPE_TRAITS_HPP
#define FWDPP_INTERNAL_TYPE_TRAITS_HPP
#include <type_traits>
#include <functional>
#include <utility>
#include <fwdpp/internal/void_t.hpp>

namespace fwdpp
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
            struct is_diploid<
                T, typename traits::internal::void_t<
                       typename T::first_type, typename T::second_type>::type>
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
                : std::integral_constant<
                      bool,
                      std::is_integral<typename T::first_type>::value
                          && std::is_same<typename T::first_type,
                                          typename T::second_type>::value>
            {
            };

            template <typename T, typename = void>
            struct is_custom_diploid : std::false_type
            {
            };

            template <typename T>
            struct is_custom_diploid<
                T, typename void_t<typename T::first_type,
                                   typename T::second_type>::type>
                : std::integral_constant<
                      bool,
                      is_diploid<T>::value
                          && !std::is_same<std::pair<typename T::first_type,
                                                     typename T::second_type>,
                                           T>::value>
            {
            };

            template <typename dipvector_t, typename gcont_t, typename mcont_t,
                      typename = void, typename = void, typename = void,
                      typename = void>
            struct fitness_fxn
            {
                using type = void;
            };

            template <typename dipvector_t, typename gcont_t, typename mcont_t>
            struct fitness_fxn<
                dipvector_t, gcont_t, mcont_t,
                typename void_t<typename dipvector_t::value_type,
                                typename gcont_t::value_type,
                                typename mcont_t::value_type>::type,
                typename std::enable_if<
                    is_diploid<typename dipvector_t::value_type>::value>::type,
                typename std::enable_if<is_haploid_genome<
                    typename gcont_t::value_type>::value>::type,
                typename std::enable_if<
                    is_mutation<typename mcont_t::value_type>::value>::type>

            {
                using type = std::function<double(
                    const typename dipvector_t::value_type &, const gcont_t &,
                    const mcont_t &)>;
            };

            template <typename ff, typename dipvector_t, typename gcont_t,
                      typename mcont_t,
                      typename ffxn_t = typename fitness_fxn<
                          dipvector_t, gcont_t, mcont_t>::type>
            struct is_fitness_fxn
                : std::integral_constant<
                      bool, !std::is_void<ffxn_t>::value
                                && std::is_convertible<ff, ffxn_t>::value>
            {
            };

            template <typename mmodel_t, typename mcont_t, typename gcont_t,
                      typename = void, typename = void, typename = void>
            struct is_mutation_model : std::false_type
            {
            };

            template <typename mmodel_t, typename mcont_t, typename gcont_t>
            struct is_mutation_model<
                mmodel_t, mcont_t, gcont_t,
                typename void_t<typename std::result_of<mmodel_t(
                    flagged_mutation_queue &, mcont_t &)>::type>::type,
                typename std::enable_if<
                    is_mutation<typename mcont_t::value_type>::value>::type,
                typename std::enable_if<is_haploid_genome<
                    typename gcont_t::value_type>::value>::type>
                : std::true_type
            {
            };

            template <typename mmodel_t, typename mcont_t, typename gcont_t>
            struct is_mutation_model<
                mmodel_t, mcont_t, gcont_t,
                typename void_t<typename std::result_of<mmodel_t(
                    flagged_mutation_queue &, typename gcont_t::value_type &,
                    mcont_t &)>::type>::type,
                typename std::enable_if<
                    is_mutation<typename mcont_t::value_type>::value>::type,
                typename std::enable_if<is_haploid_genome<
                    typename gcont_t::value_type>::value>::type>
                : std::true_type
            {
            };

            template <typename recmodel_t, typename diploid_t,
                      typename haploid_genome_t, typename mcont_t,
                      typename = void, typename = void, typename = void,
                      typename = void>
            struct is_rec_model : std::false_type
            {
            };

            template <typename recmodel_t, typename diploid_t,
                      typename haploid_genome_t, typename mcont_t>
            struct is_rec_model<
                recmodel_t, diploid_t, haploid_genome_t, mcont_t,
                typename void_t<
                    typename std::result_of<recmodel_t()>::type>::type,
                typename std::enable_if<is_diploid<diploid_t>::value>::type,
                typename std::enable_if<
                    is_haploid_genome<haploid_genome_t>::value>::type,
                typename std::enable_if<
                    is_mutation<typename mcont_t::value_type>::value>::type>
                : std::true_type
            {
            };

            template <typename recmodel_t, typename diploid_t,
                      typename haploid_genome_t, typename mcont_t>
            struct is_rec_model<
                recmodel_t, diploid_t, haploid_genome_t, mcont_t,
                typename void_t<typename std::result_of<recmodel_t(
                    const haploid_genome_t &, const haploid_genome_t &,
                    const mcont_t &)>::type>::type,
                typename std::enable_if<is_diploid<diploid_t>::value>::type,
                typename std::enable_if<
                    is_haploid_genome<haploid_genome_t>::value>::type,
                typename std::enable_if<
                    is_mutation<typename mcont_t::value_type>::value>::type>
                : std::true_type
            {
            };

            template <typename recmodel_t, typename diploid_t,
                      typename haploid_genome_t, typename mcont_t>
            struct is_rec_model<
                recmodel_t, diploid_t, haploid_genome_t, mcont_t,
                typename void_t<typename std::result_of<recmodel_t(
                    const diploid_t &, const haploid_genome_t &,
                    const haploid_genome_t &, const mcont_t &)>::type>::type,
                typename std::enable_if<is_diploid<diploid_t>::value>::type,
                typename std::enable_if<
                    is_haploid_genome<haploid_genome_t>::value>::type,
                typename std::enable_if<
                    is_mutation<typename mcont_t::value_type>::value>::type>
                : std::true_type
            {
            };

            template <typename mcont_t, typename = void> struct mutation_model
            {
                using type = void;
            };

            template <typename mcont_t>
            struct mutation_model<
                mcont_t, typename std::enable_if<is_mutation<
                             typename mcont_t::value_type>::value>::type>
            {
                using type = std::function<std::size_t(
                    flagged_mutation_queue &, mcont_t &)>;
            };

            template <typename mcont_t, typename gcont_t, typename = void,
                      typename = void>
            struct mutation_model_haploid_genome
            {
                using type = void;
            };

            template <typename mcont_t, typename gcont_t>
            struct mutation_model_haploid_genome<
                mcont_t, gcont_t,
                typename std::enable_if<
                    is_mutation<typename mcont_t::value_type>::value>::type,
                typename std::enable_if<is_haploid_genome<
                    typename gcont_t::value_type>::value>::type>
            {
                using type = std::function<std::size_t(
                    flagged_mutation_queue &,
                    const typename gcont_t::value_type &, mcont_t &)>;
            };

            template <typename diploid_t, typename mcont_t, typename gcont_t,
                      typename = void, typename = void, typename = void>
            struct mutation_model_diploid
            {
                using type = void;
            };

            template <typename diploid_t, typename mcont_t, typename gcont_t>
            struct mutation_model_diploid<
                diploid_t, mcont_t, gcont_t,
                typename std::enable_if<is_diploid<diploid_t>::value>::type,
                typename std::enable_if<
                    is_mutation<typename mcont_t::value_type>::value>::type,
                typename std::enable_if<is_haploid_genome<
                    typename gcont_t::value_type>::value>::type>
            {
                using type = std::function<std::size_t(
                    flagged_mutation_queue &, const diploid_t &,
                    const typename gcont_t::value_type &, mcont_t &)>;
            };
        } // namespace internal
    }     // namespace traits
} // namespace fwdpp

#endif
