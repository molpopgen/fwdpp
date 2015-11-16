#ifndef __FWDPP_INTERNAL_DIPLOID_FITNESS_DISPATCH_HPP__
#define __FWDPP_INTERNAL_DIPLOID_FITNESS_DISPATCH_HPP__

#include <type_traits>
#include <fwdpp/tags/diploid_tags.hpp>

namespace KTfwd {
  namespace fwdpp_internal{
    //! Standard diploids
    template<typename fitness_policy_type,
	     typename diploid_itr_t>
    inline floating_t diploid_fitness_dispatch( const fitness_policy_type & fp, const diploid_itr_t & d, std::false_type ) {
      return fp(d->first,d->second);
    }
    
    //! Custom diploids
    template<typename fitness_policy_type,
	     typename diploid_itr_t>
    inline floating_t  diploid_fitness_dispatch( const fitness_policy_type & fp, const diploid_itr_t & d, std::true_type ) { 
      return fp(d);
    }

    //! Custom diploids
    template<typename fitness_policy_type,
	     typename diploid_itr_t>
    inline floating_t  diploid_fitness_dispatch( const fitness_policy_type & fp, const diploid_itr_t & d, std::true_type,
					     typename std::enable_if< traits::is_custom_diploid_t<typename diploid_itr_t::value_type>::type >::type()) { 
      return fp(*d);
    }
    
  }
}


#endif
