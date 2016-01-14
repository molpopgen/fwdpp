#ifndef FWDPP_POLICIES_HPP
#define FWDPP_POLICIES_HPP

#include <type_traits>
#include <functional>
#include <fwdpp/type_traits.hpp>

namespace KTfwd
{
  template<typename dipvector_t,
	   typename gcont_t,
	   typename mcont_t>
  struct fitness_fxn_type
  {
    using type = typename std::conditional< traits::is_custom_diploid_t<typename dipvector_t::value_type>::value,
					    std::function<double(const typename dipvector_t::value_type &,
								 const gcont_t &,
								 const mcont_t &)>,
					    std::function<double(const typename gcont_t::value_type &,
								 const typename gcont_t::value_type &,
								 const mcont_t &)>
					    >::type;
  };

  template<typename dipvector_t,
	   typename gcont_t,
	   typename mcont_t>
  using fitness_fxn_t = typename fitness_fxn_type<dipvector_t,gcont_t,mcont_t>::type;

  //! Gives the recombination model function signature corresponding to gcont_t,mcont_t
  template<typename gcont_t,typename mcont_t>
  using recmodel_t = std::function<std::vector<double>(const typename gcont_t::value_type &,
						       const typename gcont_t::value_type &,
						       const mcont_t &)>;
}

#endif
