#ifndef FWDPP_SUGAR_DEMOGRAPHY_HPP
#define FWDPP_SUGAR_DEMOGRAPHY_HPP

#include <type_traits>
#include <fwdpp/demography.hpp>
#include <fwdpp/sugar/metapop.hpp>

namespace KTfwd
{
  template<typename mpoptype>
  void update_Ns(mpoptype & mpop)
  {
    static_assert( std::is_same<typename mpoptype::popmodel_t,
		   sugar::METAPOP_TAG>::value,
		   "mpoptype must be an object of type KTfwd::metapop or KTfwd::sugar::metapop" );
    mpop.Ns.clear();
    for( const auto & dip : mpop.diploids )
      {
	mpop.Ns.push_back(dip.size());
      }
  }
  
  template<typename mpoptype>
  int copy_pop(mpoptype & mpop,
	       const size_t i)
  {
    static_assert( std::is_same<typename mpoptype::popmodel_t,
		   sugar::METAPOP_TAG>::value,
		   "mpoptype must be an object of type KTfwd::metapop or KTfwd::sugar::metapop" );
    auto rv = copy_deme(mpop.mutations,mpop.mcounts,mpop.gametes,
			mpop.diploids,i);
    if(rv) return rv; //there was an error
    update_Ns(mpop); //update deme sizes
    return rv;
  }

  template<typename mpoptype>
  int merge_pops(mpoptype & mpop,
		 const size_t i,
		 const size_t j)
  {
    static_assert( std::is_same<typename mpoptype::popmodel_t,
		   sugar::METAPOP_TAG>::value,
		   "mpoptype must be an object of type KTfwd::metapop or KTfwd::sugar::metapop" );
    auto rv = merge_demes(mpop.diploids,i,j);
    if(rv) return rv; //there was an error
    update_Ns(mpop); //update deme sizes
    return rv;
  }

  template<typename mpoptype>
  int remove_pop(mpoptype & mpop,
		 const size_t i)
  {
    static_assert( std::is_same<typename mpoptype::popmodel_t,
		   sugar::METAPOP_TAG>::value,
		   "mpoptype must be an object of type KTfwd::metapop or KTfwd::sugar::metapop" );
    auto rv = remove_deme(mpop.mutations,mpop.mcounts,mpop.gametes,mpop.diploids,i);
    if(rv) return rv; //there was an error
    update_Ns(mpop); //update deme sizes
    return rv;    
  }

  template<typename mpoptype>
  int swap_pops(mpoptype & mpop,
		const size_t i,
		const size_t j)
  {
    static_assert( std::is_same<typename mpoptype::popmodel_t,
		   sugar::METAPOP_TAG>::value,
		   "mpoptype must be an object of type KTfwd::metapop or KTfwd::sugar::metapop" );
    return swap_demes(mpop.diploids,i,j);
  }

  template<typename mpoptype>
  int split_pop(gsl_rng * r,
		mpoptype & mpop,
		const size_t i,
		const uint_t N_new,
		const bool replacement = false)
  {
    static_assert( std::is_same<typename mpoptype::popmodel_t,
		   sugar::METAPOP_TAG>::value,
		   "mpoptype must be an object of type KTfwd::metapop or KTfwd::sugar::metapop" );
    assert(i<mpop.diploids.size());
    auto rv = split_deme(r,mpop.mutations,mpop.mcounts,mpop.gametes,mpop.diploids,i,N_new,replacement);
    if(rv) return rv;
    update_Ns(mpop);
    return rv;
  }

  template<typename mpoptype>
  int admix_pops(gsl_rng * r,
		 mpoptype & mpop,
		 const size_t i,
		 const size_t j,
		 const double pi,
		 const uint_t N_new,
		 const bool replacement = false)
  {
    static_assert( std::is_same<typename mpoptype::popmodel_t,
		   sugar::METAPOP_TAG>::value,
		   "mpoptype must be an object of type KTfwd::metapop or KTfwd::sugar::metapop" );
    auto rv = admix_demes(r,mpop.mutations,mpop.mcounts,mpop.gametes,mpop.diploids,i,j,pi,N_new,replacement);
    if(rv) return rv;
    update_Ns(mpop);
    return rv;
  }
}

#endif


