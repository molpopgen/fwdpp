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

            template <typename T, typename = void> struct is_diploid : std::false_type
            /* Fallback type.  Evaluates to std::false_type */
            {
            };

            template <typename T>
            struct is_diploid<
                T, typename traits::internal::void_t<typename T::first_type,
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
                : std::integral_constant<
                      bool, std::is_integral<typename T::first_type>::value
                                && std::is_same<typename T::first_type,
                                                typename T::second_type>::value>
            {
            };

            template <typename T, typename = void>
            struct is_custom_diploid : std::false_type
            {
            };

            template <typename T>
            struct is_custom_diploid<T, typename void_t<typename T::first_type,
                                                        typename T::second_type>::type>
                : std::integral_constant<
                      bool, is_diploid<T>::value
                                && !std::is_same<std::pair<typename T::first_type,
                                                           typename T::second_type>,
                                                 T>::value>
            {
            };

            template <typename DiploidContainerType, typename GenomeContainerType,
                      typename MutationContainerType, typename = void, typename = void,
                      typename = void, typename = void>
            struct fitness_fxn
            {
                using type = void;
            };

            template <typename DiploidContainerType, typename GenomeContainerType,
                      typename MutationContainerType>
            struct fitness_fxn<
                DiploidContainerType, GenomeContainerType, MutationContainerType,
                typename void_t<typename DiploidContainerType::value_type,
                                typename GenomeContainerType::value_type,
                                typename MutationContainerType::value_type>::type,
                typename std::enable_if<
                    is_diploid<typename DiploidContainerType::value_type>::value>::type,
                typename std::enable_if<is_haploid_genome<
                    typename GenomeContainerType::value_type>::value>::type,
                typename std::enable_if<is_mutation<
                    typename MutationContainerType::value_type>::value>::type>

            {
                using type = std::function<double(
                    const typename DiploidContainerType::value_type &,
                    const GenomeContainerType &, const MutationContainerType &)>;
            };

            template <typename FitnessFunctionType, typename DiploidContainerType,
                      typename GenomeContainerType, typename MutationContainerType,
                      typename ffxn_t =
                          typename fitness_fxn<DiploidContainerType, GenomeContainerType,
                                               MutationContainerType>::type>
            struct is_fitness_fxn
                : std::integral_constant<
                      bool,
                      !std::is_void<ffxn_t>::value
                          && std::is_convertible<FitnessFunctionType, ffxn_t>::value>
            {
            };

            template <typename MutationModel, typename MutationContainerType,
                      typename GenomeContainerType, typename = void, typename = void,
                      typename = void>
            struct is_mutation_model : std::false_type
            {
            };

            template <typename MutationModel, typename MutationContainerType,
                      typename GenomeContainerType>
            struct is_mutation_model<
                MutationModel, MutationContainerType, GenomeContainerType,
                typename void_t<typename std::result_of<MutationModel(
                    flagged_mutation_queue &, MutationContainerType &)>::type>::type,
                typename std::enable_if<is_mutation<
                    typename MutationContainerType::value_type>::value>::type,
                typename std::enable_if<is_haploid_genome<
                    typename GenomeContainerType::value_type>::value>::type>
                : std::true_type
            {
            };

            template <typename MutationModel, typename MutationContainerType,
                      typename GenomeContainerType>
            struct is_mutation_model<
                MutationModel, MutationContainerType, GenomeContainerType,
                typename void_t<typename std::result_of<MutationModel(
                    typename GenomeContainerType::value_type &, flagged_mutation_queue &,
                    MutationContainerType &)>::type>::type,
                typename std::enable_if<is_mutation<
                    typename MutationContainerType::value_type>::value>::type,
                typename std::enable_if<is_haploid_genome<
                    typename GenomeContainerType::value_type>::value>::type>
                : std::true_type
            {
            };

            template <typename RecombinationModelType, typename DiploidType,
                      typename HaploidGenomeType, typename MutationContainerType,
                      typename = void, typename = void, typename = void, typename = void>
            struct is_rec_model : std::false_type
            {
            };

            template <typename RecombinationModelType, typename DiploidType,
                      typename HaploidGenomeType, typename MutationContainerType>
            struct is_rec_model<
                RecombinationModelType, DiploidType, HaploidGenomeType,
                MutationContainerType,
                typename void_t<
                    typename std::result_of<RecombinationModelType()>::type>::type,
                typename std::enable_if<is_diploid<DiploidType>::value>::type,
                typename std::enable_if<
                    is_haploid_genome<HaploidGenomeType>::value>::type,
                typename std::enable_if<is_mutation<
                    typename MutationContainerType::value_type>::value>::type>
                : std::true_type
            {
            };

            template <typename RecombinationModelType, typename DiploidType,
                      typename HaploidGenomeType, typename MutationContainerType>
            struct is_rec_model<
                RecombinationModelType, DiploidType, HaploidGenomeType,
                MutationContainerType,
                typename void_t<typename std::result_of<RecombinationModelType(
                    const HaploidGenomeType &, const HaploidGenomeType &,
                    const MutationContainerType &)>::type>::type,
                typename std::enable_if<is_diploid<DiploidType>::value>::type,
                typename std::enable_if<
                    is_haploid_genome<HaploidGenomeType>::value>::type,
                typename std::enable_if<is_mutation<
                    typename MutationContainerType::value_type>::value>::type>
                : std::true_type
            {
            };

            template <typename RecombinationModelType, typename DiploidType,
                      typename HaploidGenomeType, typename MutationContainerType>
            struct is_rec_model<
                RecombinationModelType, DiploidType, HaploidGenomeType,
                MutationContainerType,
                typename void_t<typename std::result_of<RecombinationModelType(
                    const DiploidType &, const HaploidGenomeType &,
                    const HaploidGenomeType &, const MutationContainerType &)>::type>::
                    type,
                typename std::enable_if<is_diploid<DiploidType>::value>::type,
                typename std::enable_if<
                    is_haploid_genome<HaploidGenomeType>::value>::type,
                typename std::enable_if<is_mutation<
                    typename MutationContainerType::value_type>::value>::type>
                : std::true_type
            {
            };

            template <typename MutationContainerType, typename = void>
            struct mutation_model
            {
                using type = void;
            };

            template <typename MutationContainerType>
            struct mutation_model<
                MutationContainerType,
                typename std::enable_if<is_mutation<
                    typename MutationContainerType::value_type>::value>::type>
            {
                using type = std::function<std::size_t(flagged_mutation_queue &,
                                                       MutationContainerType &)>;
            };

            template <typename MutationContainerType, typename GenomeContainerType,
                      typename = void, typename = void>
            struct mutation_model_haploid_genome
            {
                using type = void;
            };

            template <typename MutationContainerType, typename GenomeContainerType>
            struct mutation_model_haploid_genome<
                MutationContainerType, GenomeContainerType,
                typename std::enable_if<is_mutation<
                    typename MutationContainerType::value_type>::value>::type,
                typename std::enable_if<is_haploid_genome<
                    typename GenomeContainerType::value_type>::value>::type>
            {
                using type = std::function<std::size_t(
                    const typename GenomeContainerType::value_type &,
                    flagged_mutation_queue &,
                    MutationContainerType &)>;
            };

            template <typename DiploidType, typename MutationContainerType,
                      typename GenomeContainerType, typename = void, typename = void,
                      typename = void>
            struct mutation_model_diploid
            {
                using type = void;
            };

            template <typename DiploidType, typename MutationContainerType,
                      typename GenomeContainerType>
            struct mutation_model_diploid<
                DiploidType, MutationContainerType, GenomeContainerType,
                typename std::enable_if<is_diploid<DiploidType>::value>::type,
                typename std::enable_if<is_mutation<
                    typename MutationContainerType::value_type>::value>::type,
                typename std::enable_if<is_haploid_genome<
                    typename GenomeContainerType::value_type>::value>::type>
            {
                using type = std::function<std::size_t(
                    flagged_mutation_queue &, const DiploidType &,
                    const typename GenomeContainerType::value_type &,
                    MutationContainerType &)>;
            };
        } // namespace internal
    }     // namespace traits
} // namespace fwdpp

#endif
