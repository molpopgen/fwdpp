#ifndef __FWDPP_SUGAR_MUTATION_POPGENMUT_HPP__
#define __FWDPP_SUGAR_MUTATION_POPGENMUT_HPP__

#include <fwdpp/forward_types.hpp>
#include <fwdpp/tags/tags.hpp>
#include <limits>
namespace KTfwd
{
  /*!
    \brief Mutations with selection, dominance, and tracking age of origin
    The "standard" mutation type for population genetic simulation.
    A mutation has its own selection and dominance coefficients.

    \ingroup sugar
   */
  struct popgenmut : public mutation_base
  {
    //! The generation when the mutation arose
    uint_t g;
    //! Selection coefficient
    double s;
    //! Dominance of the mutation
    double h;
    /*!
      \brief Constructor
      \param __pos Mutation position
      \param __s Selection coefficient
      \param __h Dominance coefficient
      \param __g Generation when mutation arose
      \param __n Number of copies of mutation in population
    */
    popgenmut(const double & __pos, const double & __s, const double & __h,
	      const unsigned & __g) noexcept
      : mutation_base(__pos,(__s==0.) ? true : false),g(__g),s(__s),h(__h)
    {
    }
    popgenmut( tags::extinct ) noexcept
    : mutation_base(std::numeric_limits<double>::quiet_NaN(),true),
      g(std::numeric_limits<unsigned>::max()),
      s(std::numeric_limits<double>::quiet_NaN()),
      h(std::numeric_limits<double>::quiet_NaN())
    {
    }
  };
}
#endif
