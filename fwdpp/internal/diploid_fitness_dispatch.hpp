#ifndef __FWDPP_INTERNAL_DIPLOID_FITNESS_DISPATCH_HPP__
#define __FWDPP_INTERNAL_DIPLOID_FITNESS_DISPATCH_HPP__

#include <type_traits>
#include <fwdpp/tags/diploid_tags.hpp>

namespace KTfwd {
  namespace fwdpp_internal{
    template<typename fitness_policy_type,
	     typename diploid_itr_t>
    double diploid_fitness_dispatch( const fitness_policy_type & fp, const diploid_itr_t & d, const KTfwd::tags::standard_diploid_t & ) {
      return fp(d->first,d->second);
    }

    template<typename fitness_policy_type,
	     typename diploid_itr_t>
    double  diploid_fitness_dispatch( const fitness_policy_type & fp, const diploid_itr_t & d, const KTfwd::tags::custom_diploid_t & ) { 
      return fp(d);
    }
    
  }
}


#endif
