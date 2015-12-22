#ifndef __FWDPP_TYPE_TRAITS_HPP__
#define __FWDPP_TYPE_TRAITS_HPP__

#include <type_traits>
#include <functional>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/tags/diploid_tags.hpp>
#include <fwdpp/internal/recycling.hpp>
#include <fwdpp/internal/gamete_lookup_table.hpp>
#include <fwdpp/internal/mutation_internal.hpp>

namespace KTfwd {
  namespace traits {
    //! Evaluates to std::true_type if T inherits from KTfwd::mutation_base
    template<typename T>
    using is_mutation_t = typename std::is_base_of<KTfwd::mutation_base,T>::type;

    //! Evaluates to std::true_type if T publicly inherits from KTfwd::tags::custom_diploid_t
    template<typename T>
    using is_custom_diploid_t = typename std::is_base_of<KTfwd::tags::custom_diploid_t,T>::type;

    template<typename T>
    struct has_gamete_tag
    {
    private:
      typedef char                      yes;
      typedef struct { char array[2]; } no;
      
      template<typename C> static yes test(typename C::gamete_tag*);
      template<typename C> static no  test(...);
    public:
      static const bool value = sizeof(test<T>(0)) == sizeof(yes);
    };

    template<class T>
    struct void_t {
      typedef void type;
    };

    template<class T, class U = void>
    struct has_gamete_tag2 : std::false_type {
    };

    template<class T>
    struct has_gamete_tag2<T, typename void_t<typename T::gamete_tag>::type > : std::true_type
    {
    };
    
    //! Gives the "recycling bin" type corresponding to list_t
    template<typename list_t>
    struct recycling_bin_t
    {
      using type = KTfwd::fwdpp_internal::recycling_bin_t<typename list_t::iterator>;
    };

    //! Gives the "gamete lookup table" type corresponding to list_t
    template<typename list_t>
    struct gamete_lookup_t
    {
      static_assert( has_gamete_tag<typename list_t::value_type>::value , "foo" );
      static_assert( has_gamete_tag2<typename list_t::value_type>::value , "foo" );
      using type = typename std::result_of<decltype(&fwdpp_internal::gamete_lookup_table<list_t>)(list_t *)>::type;
    };

    /*! 
      Gives the mutation model function signature corresponding to mlist_t.

      Applies to mutation policies that only take recycling bins and  mlist_t *
      as arguments
    */
    template<typename mlist_t>
    struct mmodel_t
    {
      static_assert( typename is_mutation_t<typename mlist_t::value_type>::type(),
		     "mlist_t::value_type must be derived from KTfwd::mutation_base" );
      using type = std::function<typename mlist_t::iterator(typename recycling_bin_t<mlist_t>::type &,mlist_t *)>;
    };

    /*!
      Gives mutation model function signature for models requiring gametes as arguments
    */
    template<typename mlist_t,typename glist_t>
    struct mmodel_gamete_t
    {
      using type = std::function<typename mlist_t::iterator(typename recycling_bin_t<mlist_t>::type &,
							    typename glist_t::value_type &,
							    mlist_t *)>;
    };

    //! Check that a mutation model type is valid.
    template<typename mmodel_t,typename mlist_t,typename glist_t>
    struct valid_mutation_model
    {
      using queue_t = typename recycling_bin_t<mlist_t>::type;
      using gamete_t = typename glist_t::value_type;
      //http://stackoverflow.com/questions/11470802/stdresult-of-simple-function
      //http://stackoverflow.com/questions/2763824/decltype-result-of-or-typeof
      using result_type = typename std::result_of<decltype(&fwdpp_internal::mmodel_dispatcher<mmodel_t,gamete_t,mlist_t,queue_t>)(mmodel_t &,gamete_t &,mlist_t *,queue_t &)>::type;
      using type = typename std::is_same<result_type,typename mlist_t::iterator>::type;
    };
      
    //! Gives the recombination model function signature corresponding to glist_t
    template<typename glist_t>
    struct recmodel_t
    {
      using type = std::function<unsigned(typename glist_t::iterator &,typename glist_t::iterator &,
					  typename gamete_lookup_t<glist_t>::type &,
					  typename recycling_bin_t<glist_t>::type &) >;
    };

  }
} 
#endif
