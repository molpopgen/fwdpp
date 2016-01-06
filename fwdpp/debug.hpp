#ifndef __KTFWD_DEBUG_HPP__
#define __KTFWD_DEBUG_HPP__

#ifndef NDEBUG
#include <algorithm>
#include <numeric>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/type_traits.hpp>

namespace KTfwd
{

  /*! \brief Returns true if the sum of counts in gametes equals twoN, false otherwise
    Returns true if the sum of counts in gametes equals twoN, false otherwise
   */
  template<typename gcont_t>
  bool check_sum(const gcont_t & gametes, const unsigned & twoN)
  {
    static_assert( typename traits::is_gamete_t<typename gcont_t::value_type>::type(),
		   "gcont_t::value_type must be a valid gamete type" );
    return ( std::accumulate( gametes.cbegin(),
			      gametes.cend(),0u,
			      [](unsigned & __u,
				 const typename gcont_t::value_type & __g) { 
				return __u + __g.n; 
			      } ) == twoN );
  }

  /*! \brief Returns true if the sum of counts in gametes equals twoN, false otherwise
    Returns true if the sum of counts in gametes equals twoN, false otherwise
   */
  template<typename gcont_t>
  bool check_sum(const gcont_t * gametes, const unsigned & twoN)
  {
    return check_sum(*gametes,twoN);
  }

  template<typename gamete_t,
	   typename mcont_t>
  bool gamete_is_sorted_n( const gamete_t & g,
			   const mcont_t & m )
  {
    return std::is_sorted(g.mutations.begin(),
			  g.mutations.end(),
			  [&m](const size_t i, const size_t j)
			  {
			    return m[i].pos <= m[j].pos;
			  });
  }

  template<typename gamete_t,
	   typename mcont_t>
  bool gamete_is_sorted_s( const gamete_t & g,
			   const mcont_t & m )
  {
    return std::is_sorted(g.smutations.begin(),
			  g.smutations.end(),
			  [&m](const size_t i, const size_t j)
			  {
			    return m[i].pos <= m[j].pos;
			  });
  }

}

#endif
#endif
