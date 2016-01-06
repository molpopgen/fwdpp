#ifndef __FWDPP_INTERNAL_DIPLOID_FITNESS_DISPATCH_HPP__
#define __FWDPP_INTERNAL_DIPLOID_FITNESS_DISPATCH_HPP__

#include <type_traits>
#include <fwdpp/tags/diploid_tags.hpp>

namespace KTfwd {
  namespace fwdpp_internal{
    //! Standard diploids
    template<typename fitness_policy_type,
	     typename diploid_itr_t,
	     typename gcont_t,
	     typename mcont_t>
    inline double diploid_fitness_dispatch( const fitness_policy_type & fp, const diploid_itr_t & d,
					    const gcont_t & gametes,const mcont_t & mutations,
					    std::false_type ) {
      return fp(gametes[d.first],gametes[d.second],mutations);
    }
    
    //! Custom diploids
    template<typename fitness_policy_type,
	     typename diploid_itr_t,
	     typename gcont_t,
	     typename mcont_t>
    inline double  diploid_fitness_dispatch( const fitness_policy_type & fp, const diploid_itr_t & d,
					     const gcont_t & gametes,const mcont_t & mutations,
					     std::true_type ) { 
      return fp(d,gametes,mutations);
    }

    //! Custom diploids
    template<typename fitness_policy_type,
	     typename diploid_itr_t,
	     typename gcont_t,
	     typename mcont_t>
    inline double  diploid_fitness_dispatch( const fitness_policy_type & fp, const diploid_itr_t & d,
					     const gcont_t & gametes,const mcont_t & mutations,
					     std::true_type,
					     typename std::enable_if< traits::is_custom_diploid_t<typename diploid_itr_t::value_type>::type >::type()) { 
      return fp(d,gametes,mutations);
    }
    
  }
}


#endif
