#ifndef FWDPP_SUGAR_GENERALMUT_HPP
#define FWDPP_SUGAR_GENERALMUT_HPP

#include <array>
#include <memory>

#include <fwdpp/forward_types.hpp>

namespace KTfwd
{
  /*!
    \brief Mutation type allowing arbitray number of "s,h" pairs representing effect sizes,dominance.
    \ingroup sugar
  */
  template<std::size_t N>
  struct generalmut : public KTfwd::mutation_base
  {
    using array_t = std::array<double,N>;
    //! Selection coefficients and/or effect sizes
    array_t s;
    //! Dominances associated w/values in 's'.
    array_t h;
    //! Generation when mutation arose
    unsigned g;
    //! Constructor
    generalmut( array_t __s,
		array_t __h,
		double pos,unsigned n,unsigned gen ) : KTfwd::mutation_base(std::move(pos),std::move(n),
									    (std::find_if(std::begin(__s),std::end(__s),[](const double d) { return d != 0.; }) == std::end(__s))),
						       s(std::move(__s)),h(std::move(__h)),g(std::move(gen))
    {
    }
  };
}

#endif
