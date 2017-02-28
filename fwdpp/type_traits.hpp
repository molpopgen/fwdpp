#ifndef __FWDPP_TYPE_TRAITS_HPP__
#define __FWDPP_TYPE_TRAITS_HPP__

#include <type_traits>
#include <functional>
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

        //! Wraps a static constant allowing a test that T is a mutation
        template <typename T>
        struct is_mutation
            : std::integral_constant<bool,
                                     std::is_base_of<KTfwd::mutation_base,
                                                     T>::value>
        {
        };

		//! Gives the "recycling bin" type corresponding to cont_t
        template <typename cont_t> struct recycling_bin_type
        {
            using type = KTfwd::fwdpp_internal::recycling_bin_t<
                typename cont_t::size_type>;
        };

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

        template <typename T>
        using is_multilocus_diploid
            = traits::internal::is_multilocus_diploid<T>;

        //! Wraps a static constant allowing a test that T is a custom diploid
        template <typename T>
        using is_custom_diploid = traits::internal::is_custom_diploid<T>;

        //! Evaluates to KTfwd::traits::is_diploid<T>::type
        template <typename T>
        using is_diploid_t = typename is_diploid<T>::type;

        //! Evaluates to KTfwd::traits::is_custom_diploid<T>::type
        template <typename T>
        using is_custom_diploid_t = typename is_custom_diploid<T>::type;

#if __cplusplus >= 201402L
        template <typename T>
        constexpr bool is_diploid_v = is_diploid<T>::value;

        template <typename T>
        constexpr bool is_custom_diploid_v = is_custom_diploid<T>::value;
#endif

        //! Wraps a static constant to test that mmodel_t is a valid mutation
        //! model/policy
        template <typename mmodel_t, typename mcont_t, typename gcont_t>
        using is_mutation_model
            = traits::internal::is_mutation_model<mmodel_t, mcont_t, gcont_t>;

        //! Wraps a static constant to test that recmodel_t is a valid mutation
        //! model/policy
        template <typename recmodel_t, typename gamete_t, typename mcont_t>
        struct is_rec_model
            : std::
                  integral_constant<bool,
                                    std::is_same<
                                        typename std::result_of<recmodel_t(
                                            const gamete_t &, const gamete_t &,
                                            const mcont_t &)>::type,
                                        std::vector<double>>::value>
        {
        };

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

        //! Evaulates to fitness_fxn<dipvector_t,gcont_t,mcont_t>::type
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

        //! Gives the recombination model function signature corresponding
        //! to
        //! gcont_t,mcont_t
        template <typename gcont_t, typename mcont_t>
        using recmodel_t = std::function<std::vector<double>(
            const typename gcont_t::value_type &,
            const typename gcont_t::value_type &, const mcont_t &)>;

        /*!
          Gives the mutation model function signature corresponding to mcont_t.

          Applies to mutation policies that only take recycling bins and
          mcont_t *
          as arguments
        */
        template <typename mcont_t>
        using mmodel_t = std::function<std::size_t(recycling_bin_t<mcont_t> &,
                                                   mcont_t &)>;

        /*!
          Gives mutation model function signature for models requiring gametes
          as arguments
        */
        template <typename mcont_t, typename gcont_t>
        using mmodel_gamete_t
            = std::function<std::size_t(recycling_bin_t<mcont_t> &,
                                        typename gcont_t::value_type &,
                                        mcont_t &)>;
    }
}
#endif
