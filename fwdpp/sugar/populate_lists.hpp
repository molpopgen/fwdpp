#ifndef FWDPP_SUGAR_POPULATE_LISTS_HPP
#define FWDPP_SUGAR_POPULATE_LISTS_HPP

#include <fwdpp/util.hpp>
#include <fwdpp/tags/tags.hpp>
namespace KTfwd
{
  /*!
    Add "reyclable" elements to a population object
    
    \param t A population object from fwdpp's "sugar" sub-library
    \param ngametes The final desired size of t.gametes
    \param nmutations The final desired size of t.mutations

    Recommended values for ngametes and nmutations are 2*N and 
    ceil(log(2N)*4Nu + (2/3)*4Nu), respectively.  "N" should be
    treated as the max N possible in a simulation (e.g., accounting for
    size changes).  Further, "u" should be the total mutation rate.

    The expression log(2N)*4Nu + (2/3)*4Nu is the expected number of 
    mutations in a Wright-Fisher population of size N, assuming an 
    infinitely-many sites mutation scheme.  This expression is found in
    Ewens, 2nd ed, p. 298, eqn. 9.19.  The 2/3 is really Euler's Gamma,
    but this is close enough.  This value represents a good "first guess"
    at how big a mutation list can get.

    \note Implememented via KTfwd::add_elements and KTfwd::tags::extinct
   */
  template<typename T>
  inline void add_recyclable( T & t, size_t ngametes, size_t nmutations )
  {
    add_elements(t.gametes,ngametes,tags::extinct());
    add_elements(t.mutations,nmutations,tags::extinct());
  }
}

#endif
