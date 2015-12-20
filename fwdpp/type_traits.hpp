#ifndef __FWDPP_TYPE_TRAITS_HPP__
#define __FWDPP_TYPE_TRAITS_HPP__

#include <functional>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/tags/diploid_tags.hpp>
#include <fwdpp/internal/recycling.hpp>
#include <fwdpp/internal/gamete_lookup_table.hpp>
namespace KTfwd {
  namespace traits {
    //! Evaluates to std::true_type if T inherits from KTfwd::mutation_base
    template<typename T>
    using is_mutation_t = typename std::is_base_of<KTfwd::mutation_base,T>::type;
    
    //! Evaluates to std::true_type if T publicly inherits from KTfwd::tags::custom_diploid_t
    template<typename T>
    using is_custom_diploid_t = typename std::is_base_of<KTfwd::tags::custom_diploid_t,T>::type;
    
    //! Gives the "recycling bin" type corresponding to list_t
    template<typename list_t>
    struct recycling_bin_t
    {
      using type = KTfwd::fwdpp_internal::recycling_bin_t<typename list_t::iterator>;
    };

    //! Gives the mutation model function signature corresponding to mlist_t
    template<typename mlist_t>
    struct mmodel_t
    {
      static_assert( typename is_mutation_t<typename mlist_t::value_type>::type(),
		     "mlist_t::value_type must be derived from KTfwd::mutation_base" );
      using type = std::function<typename mlist_t::iterator(typename recycling_bin_t<mlist_t>::type &,mlist_t *)>;
    };

    //! Gives the recombination model function signature corresponding to glist_t
    template<typename glist_t>
    struct recmodel_t
    {
      using type = std::function<unsigned(typename glist_t::iterator &,typename glist_t::iterator &,
					  typename KTfwd::fwdpp_internal::gamete_lookup<glist_t> &,
					  typename recycling_bin_t<glist_t>::type &) >;
    };

  }
} 
#endif
