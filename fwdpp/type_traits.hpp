#ifndef __FWDPP_TYPE_TRAITS_HPP__
#define __FWDPP_TYPE_TRAITS_HPP__

#include <type_traits>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/internal/void_t.hpp>
#include <fwdpp/internal/recycling.hpp>
#include <fwdpp/internal/mutation_internal.hpp>

namespace KTfwd
{
    namespace traits
    {
        //! Wraps a static constant allowing a test that T is a gamete
        template <typename T>
        struct is_gamete
            : std::integral_constant<bool, traits::internal::
                                               has_gamete_tag<T>::value>
        {
        };

        //! Convenience wrapper for KTfwd::traits::is_gamete<T>::type.
        template <typename T> using is_gamete_t = typename is_gamete<T>::type;

        //! Wraps a static constant allowing a test that T is a mutation
        template <typename T>
        struct is_mutation
            : std::integral_constant<bool,
                                     std::is_base_of<KTfwd::mutation_base,
                                                     T>::value>
        {
        };

        //! Convenience wrapper for KTfwd::traits::is_mutation<T>::Type.
        template <typename T>
        using is_mutation_t = typename is_mutation<T>::type;

        //! Gives the "recycling bin" type corresponding to cont_t
        template <typename cont_t> struct recycling_bin_type
        {
            using type = KTfwd::fwdpp_internal::recycling_bin_t<
                typename cont_t::size_type>;
        };

        // Evaluates to KTfwd::traits::recycling_bin_type<T>::type
        template <typename T>
        using recycling_bin_t = typename recycling_bin_type<T>::type;
    }
}

#include <fwdpp/internal/type_traits.hpp>

namespace KTfwd
{
    namespace traits
    {
        //! Wraps a static constant allowing a test that T is a diploid
        template <typename T>
        using is_diploid = traits::internal::is_diploid<T>;

        //! Wraps a static constant allowing a test that T is a multi-locus
        //! diploid.
        template <typename T>
        using is_multilocus_diploid
            = traits::internal::is_multilocus_diploid<T>;

        //! Convenience wrapper for
        //! KTfwd::traits::is_multilocus_diploid<T>::type
        template <typename T>
        using is_multilocus_diploid_t =
            typename is_multilocus_diploid<T>::type;

        //! Wraps a static constant allowing a test that T is a custom diploid
        template <typename T>
        using is_custom_diploid = traits::internal::is_custom_diploid<T>;

        //! Convenience wrapper for KTfwd::traits::is_diploid<T>::type
        template <typename T>
        using is_diploid_t = typename is_diploid<T>::type;

        //! Convenience wrapper for KTfwd::traits::is_custom_diploid<T>::type
        template <typename T>
        using is_custom_diploid_t = typename is_custom_diploid<T>::type;

        /*!
         *  Wraps a static constant to test that mmodel_t is a
         *  valid mutation policy
         */
        template <typename mmodel_t, typename mcont_t, typename gcont_t>
        using is_mutation_model
            = traits::internal::is_mutation_model<mmodel_t, mcont_t, gcont_t>;

        /*!
         * Convenience wrapper for
         * KTfwd::traits::is_mutation_model<mmodel_t,mcont_t,gcont_t>::type
         */
        template <typename mmodel_t, typename mcont_t, typename gcont_t>
        using is_mutation_model_t =
            typename is_mutation_model<mmodel_t, mcont_t, gcont_t>::type;

        /*!
         * Wraps a static constant to test that recmodel_t is a valid
         * recombination function/policy
         */
        template <typename recmodel_t, typename gamete_t, typename mcont_t>
        using is_rec_model
            = traits::internal::is_rec_model<recmodel_t, gamete_t, mcont_t>;

        /*!
         * Convenience wrapper for
         * KTfwd::traits::is_rec_model<recmodel_t,gamete_c,mcont_t>::type
         */
        template <typename recmodel_t, typename gamete_t, typename mcont_t>
        using is_rec_model_t =
            typename is_rec_model<recmodel_t, gamete_t, mcont_t>::type;

        /*!
         * Defines a struct with a single member typedef called type.
         * If type is void, then one or more of dipvector_t, gcont_t, and/or
         * mcont_t are not valid input types for a fitness function.
         * Otherwise, type will evaluate to
         * std::function<void(const dipvector_t::value_type,const gcont_t,const
         * mcont_t)>
         */
        template <typename dipvector_t, typename gcont_t, typename mcont_t>
        using fitness_fxn
            = traits::internal::fitness_fxn<dipvector_t, gcont_t, mcont_t>;

        //! Convenience wrapper for
        //! fitness_fxn<dipvector_t,gcont_t,mcont_t>::type
        template <typename dipvector_t, typename gcont_t, typename mcont_t>
        using fitness_fxn_t =
            typename fitness_fxn<dipvector_t, gcont_t, mcont_t>::type;

        /*!
         * Wrap a static constant if ff is a valid fitness function
         * whose argument types are const references to the other
         * three template parameters
         */
        template <typename ff, typename dipvector_t, typename gcont_t,
                  typename mcont_t>
        using is_fitness_fxn
            = traits::internal::is_fitness_fxn<ff, dipvector_t, gcont_t,
                                               mcont_t>;

        /*!
         * Conveneince wrapper for
         * KTfwd::traits::is_fitness_fxn<ff,dipvector_t,gcont_t,mcont_t>::type
         */
        template <typename ff, typename dipvector_t, typename gcont_t,
                  typename mcont_t>
        using is_fitness_fxn_t =
            typename is_fitness_fxn<ff, dipvector_t, gcont_t, mcont_t>::type;

        // clang-format off
        /*!
		 * Infers the signature of a recombination function compatible
		 * with the template type parameters.
		 *
		 * If such a function signature cannot be inferred, this
		 * evaluates to void.
		 */
        // clang-format on
        template <typename gcont_t_or_gamete_t, typename mcont_t>
        using recmodel_t =
            typename traits::internal::recmodel_t<gcont_t_or_gamete_t,
                                                  mcont_t>::type;

        // clang-format off
        /*!
		 * Gives the mutation model function signature corresponding to
         * mcont_t.
		 * Applies to mutation policies that only take recycling bins
         * and
		 * mcont_t as arguments.
		 *
		 * If mcont_t is not a container of mutations, then mmodel_t
		 * will evaulate to void.
		 *
		 * Otherwise, it will evaluate to
		 * std::function<std::size_t(recycling_bin_t<mcont_t> &,mcont_t
         * &)>;
         */
        // clang-format on
        template <typename mcont_t>
        using mmodel_t = typename traits::internal::mmodel_t<mcont_t>::type;

        // clang-format off
		/*!
         * Gives mutation model function signature for models requiring
         * gametes
         * as arguments. If mcont_t is not a container of mutations and/or
         * gcont_t is
		 * not a container of gametes, them mmodel_gamete_t will
         * evaluate to
		 * void.
		 *
		 * Otherwise, it will evaluate to
		 * std::function<std::size_t(recycling_bin_t<mcont_t> &,
		 * typename gcont_t::value_type &,
	     * mcont_t &)>;
        */
        // clang-format on
        template <typename mcont_t, typename gcont_t>
        using mmodel_gamete_t =
            typename traits::internal::mmodel_gamete_t<mcont_t, gcont_t>::type;

/*! \defgroup Cpp14
 * \brief C++14 features
 */

/* Template variables are provided for easier type trait
 * checking.  This feature requires compiling with the C++14
 * standard
 */
#if __cplusplus >= 201402L
        //! \ingroup Cpp14
        template <typename T>
        constexpr bool is_diploid_v = is_diploid<T>::value;

        //! \ingroup Cpp14
        template <typename T>
        constexpr bool is_custom_diploid_v = is_custom_diploid<T>::value;

        //! \ingroup Cpp14
        template <typename T>
        constexpr bool is_multilocus_diploid_v
            = is_multilocus_diploid<T>::value;

        //! \ingroup Cpp14
        template <typename T>
        constexpr bool is_mutation_v = is_mutation<T>::value;

        //! \ingroup Cpp14
        template <typename T> constexpr bool is_gamete_v = is_gamete<T>::value;

        //! \ingroup Cpp14
        template <typename mmodel_t, typename mcont_t, typename gcont_t>
        constexpr bool is_mutation_model_v
            = is_mutation_model<mmodel_t, mcont_t, gcont_t>::value;

        //! \ingroup Cpp14
        template <typename recmodel_t, typename gamete_t, typename mcont_t>
        constexpr bool is_rec_model_v
            = is_rec_model<recmodel_t, gamete_t, mcont_t>::value;

        //! \ingroup Cpp14
        template <typename ff, typename dipvector_t, typename gcont_t,
                  typename mcont_t>
        constexpr bool is_fitness_fxn_v
            = is_fitness_fxn<ff, dipvector_t, gcont_t, mcont_t>::value;
#endif
    }
}
#endif
