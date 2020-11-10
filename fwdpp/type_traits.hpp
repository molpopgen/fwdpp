#ifndef __FWDPP_TYPE_TRAITS_HPP__
#define __FWDPP_TYPE_TRAITS_HPP__

#include <type_traits>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/internal/void_t.hpp>
#include <fwdpp/simfunctions/recycling.hpp>
#include <fwdpp/internal/mutation_internal.hpp>

namespace fwdpp
{
    namespace traits
    {
        //! Wraps a static constant allowing a test that T is a haploid_genome
        template <typename T, typename = void> struct is_haploid_genome : std::false_type
        {
        };

        template <typename T>
        struct is_haploid_genome<
            T, typename traits::internal::void_t<typename T::mutation_container>::type>
            : std::integral_constant<
                  bool,
                  std::is_integral<typename T::mutation_container::value_type>::value>
        {
        };

        //! Convenience wrapper for fwdpp::traits::is_haploid_genome<T>::type.
        template <typename T>
        using is_HaploidGenomeType = typename is_haploid_genome<T>::type;

        //! Wraps a static constant allowing a test that T is a mutation
        template <typename T>
        struct is_mutation
            : std::integral_constant<bool,
                                     std::is_base_of<fwdpp::mutation_base, T>::value>
        {
        };

        //! Convenience wrapper for fwdpp::traits::is_mutation<T>::Type.
        template <typename T> using is_mutation_t = typename is_mutation<T>::type;
    } // namespace traits
} // namespace fwdpp

#include <fwdpp/internal/type_traits.hpp>

namespace fwdpp
{
    namespace traits
    {
        //! Wraps a static constant allowing a test that T is a diploid
        template <typename T> using is_diploid = traits::internal::is_diploid<T>;

        //! Wraps a static constant allowing a test that T is a custom diploid
        template <typename T>
        using is_custom_diploid = traits::internal::is_custom_diploid<T>;

        //! Convenience wrapper for fwdpp::traits::is_diploid<T>::type
        template <typename T> using is_DiploidType = typename is_diploid<T>::type;

        //! Convenience wrapper for fwdpp::traits::is_custom_diploid<T>::type
        template <typename T>
        using is_custom_diploid_t = typename is_custom_diploid<T>::type;

        /*!
         *  Wraps a static constant to test that MutationModel is a
         *  valid mutation policy
         */
        template <typename MutationModel, typename MutationContainerType,
                  typename GenomeContainerType>
        using is_mutation_model
            = traits::internal::is_mutation_model<MutationModel, MutationContainerType,
                                                  GenomeContainerType>;

        /*!
         * Convenience wrapper for
         * fwdpp::traits::is_mutation_model<MutationModel,MutationContainerType,GenomeContainerType>::type
         */
        template <typename MutationModel, typename MutationContainerType,
                  typename GenomeContainerType>
        using is_mutation_model_t =
            typename is_mutation_model<MutationModel, MutationContainerType,
                                       GenomeContainerType>::type;

        /*!
         * Wraps a static constant to test that RecombinationModelType is a valid
         * recombination function/policy
         */
        template <typename RecombinationModelType, typename DiploidType,
                  typename HaploidGenomeType, typename MutationContainerType>
        using is_rec_model
            = traits::internal::is_rec_model<RecombinationModelType, DiploidType,
                                             HaploidGenomeType, MutationContainerType>;

        /*!
         * Convenience wrapper for
         * fwdpp::traits::is_rec_model<RecombinationModelType,haploid_genome_c,MutationContainerType>::type
         */
        template <typename RecombinationModelType, typename DiploidType,
                  typename HaploidGenomeType, typename MutationContainerType>
        using is_rec_model_t =
            typename is_rec_model<RecombinationModelType, DiploidType, HaploidGenomeType,
                                  MutationContainerType>::type;

        /*!
         * Defines a struct with a single member typedef called type.
         * If type is void, then one or more of DiploidContainerType, GenomeContainerType, and/or
         * MutationContainerType are not valid input types for a fitness function.
         * Otherwise, type will evaluate to
         * std::function<void(const DiploidContainerType::value_type,const GenomeContainerType,const
         * MutationContainerType)>
         */
        template <typename DiploidContainerType, typename GenomeContainerType,
                  typename MutationContainerType>
        using fitness_fxn
            = traits::internal::fitness_fxn<DiploidContainerType, GenomeContainerType,
                                            MutationContainerType>;

        //! Convenience wrapper for
        //! fitness_fxn<DiploidContainerType,GenomeContainerType,MutationContainerType>::type
        template <typename DiploidContainerType, typename GenomeContainerType,
                  typename MutationContainerType>
        using fitness_fxn_t =
            typename fitness_fxn<DiploidContainerType, GenomeContainerType,
                                 MutationContainerType>::type;

        /*!
         * Wrap a static constant if FitnessFunctionType is a valid fitness function
         * whose argument types are const references to the other
         * three template parameters
         */
        template <typename FitnessFunctionType, typename DiploidContainerType,
                  typename GenomeContainerType, typename MutationContainerType>
        using is_fitness_fxn
            = traits::internal::is_fitness_fxn<FitnessFunctionType, DiploidContainerType,
                                               GenomeContainerType,
                                               MutationContainerType>;

        /*!
         * Conveneince wrapper for
         * fwdpp::traits::is_fitness_fxn<FitnessFunctionType,DiploidContainerType,GenomeContainerType,MutationContainerType>::type
         */
        template <typename FitnessFunctionType, typename DiploidContainerType,
                  typename GenomeContainerType, typename MutationContainerType>
        using is_fitness_fxn_t =
            typename is_fitness_fxn<FitnessFunctionType, DiploidContainerType,
                                    GenomeContainerType, MutationContainerType>::type;

        // clang-format off
        /*!
		 * Gives the mutation model function signature corresponding to
         * MutationContainerType.
		 * Applies to mutation policies that only take recycling bins
         * and
		 * MutationContainerType as arguments.
		 *
		 * If MutationContainerType is not a container of mutations, then MutationModel
		 * will evaulate to void.
		 *
		 * Otherwise, it will evaluate to
		 * std::function<std::size_t(recycling_bin_t<MutationContainerType> &,MutationContainerType
         * &)>;
         */
        // clang-format on
        template <typename MutationContainerType>
        using mutation_model =
            typename traits::internal::mutation_model<MutationContainerType>::type;

        // clang-format off
		/*!
         * Gives mutation model function signature for models requiring
         * haploid_genomes
         * as arguments. If MutationContainerType is not a container of mutations and/or
         * GenomeContainerType is
		 * not a container of haploid_genomes, them mmodel_HaploidGenomeType will
         * evaluate to
		 * void.
		 *
		 * Otherwise, it will evaluate to
		 * std::function<std::size_t(recycling_bin_t<MutationContainerType> &,
		 * const typename GenomeContainerType::value_type &,
	     * MutationContainerType &)>;
        */
        // clang-format on
        template <typename MutationContainerType, typename GenomeContainerType>
        using mutation_model_haploid_genome =
            typename traits::internal::mutation_model_haploid_genome<
                MutationContainerType, GenomeContainerType>::type;

        // clang-format off
		/*!
         * Gives mutation model function signature for models requiring
         * diploids
         * as arguments. If MutationContainerType is not a container of mutations and/or
         * GenomeContainerType is
		 * not a container of haploid_genomes, them mmodel_HaploidGenomeType will
         * evaluate to
		 * void.
		 *
		 * Otherwise, it will evaluate to
		 * std::function<std::size_t(recycling_bin_t<MutationContainerType> &,
		 * const DiploidType &,const typename GenomeContainerType::value_type &,
	     * MutationContainerType &)>;
        */
        // clang-format on
        template <typename DiploidType, typename MutationContainerType,
                  typename GenomeContainerType>
        using mutation_model_diploid = typename traits::internal::mutation_model_diploid<
            DiploidType, MutationContainerType, GenomeContainerType>::type;

        /*! \defgroup Cpp14 C++14 features
         * \brief C++14 features
         */

        /* Template variables are provided for easier type trait
         * checking.  This feature requires compiling with the C++14
         * standard.
         */
        //! \ingroup Cpp14
        template <typename T> constexpr bool is_diploid_v = is_diploid<T>::value;

        //! \ingroup Cpp14
        template <typename T>
        constexpr bool is_custom_diploid_v = is_custom_diploid<T>::value;

        //! \ingroup Cpp14
        template <typename T> constexpr bool is_mutation_v = is_mutation<T>::value;

        //! \ingroup Cpp14
        template <typename T>
        constexpr bool is_haploid_genome_v = is_haploid_genome<T>::value;

        //! \ingroup Cpp14
        template <typename MutationModel, typename MutationContainerType,
                  typename GenomeContainerType>
        constexpr bool is_mutation_model_v
            = is_mutation_model<MutationModel, MutationContainerType,
                                GenomeContainerType>::value;

        //! \ingroup Cpp14
        template <typename RecombinationModelType, typename DiploidType,
                  typename HaploidGenomeType, typename MutationContainerType>
        constexpr bool is_rec_model_v
            = is_rec_model<RecombinationModelType, DiploidType, HaploidGenomeType,
                           MutationContainerType>::value;

        //! \ingroup Cpp14
        template <typename ff, typename DiploidContainerType,
                  typename GenomeContainerType, typename MutationContainerType>
        constexpr bool is_fitness_fxn_v
            = is_fitness_fxn<ff, DiploidContainerType, GenomeContainerType,
                             MutationContainerType>::value;
    } // namespace traits
} // namespace fwdpp
#endif
