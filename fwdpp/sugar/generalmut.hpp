#ifndef FWDPP_SUGAR_GENERALMUT_HPP
#define FWDPP_SUGAR_GENERALMUT_HPP

#include <array>
#include <algorithm>
#include <memory>
#include <limits>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/tags/tags.hpp>

namespace KTfwd
{
    /*!
      \brief Mutation type allowing arbitray number of "s,h" pairs representing
      effect sizes,dominance.
      \ingroup sugar

      This type represents a mutation with a number of 's' and 'h' parameters
      that are fixed
      at compile time.

      For example, generalmut<2> could represent separate effects on fitness
      and a trait, while
      generalmut<4> could represent separate effects on fitness and on a trait
      in males and females.

      Other example use cases could be generalmut<N> representing
      fitness/dominance value pairs in
      different 'patches' of an environment.

      The advanatage of this type is that std::array<double,N> is a very fast
      type (often faster
      than a std::vector<double> of the same size).  Further, this type can be
      filled via compile-time
      meta-programming techniques.

      The disadvantage of this type is that it is not easy to have the number
      of 's,h' pairs be flexible
      at run-time.  For that use case, see KTfwd::generalmut_vec, which uses
      std::vector<double> instead of
      std::array to store the 's,h' values.
    */
    template <std::size_t N> struct generalmut : public mutation_base
    {
        using array_t = std::array<double, N>;
        //! Selection coefficients and/or effect sizes
        array_t s;
        //! Dominances associated w/values in 's'.
        array_t h;
        //! Generation when mutation arose
        uint_t g;
        //! Constructor
        generalmut(array_t __s, array_t __h, double pos, uint_t gen)
            : mutation_base(
                  std::move(pos),
                  // Mutation is neutral i.f.f. all values in __s == 0.
                  (std::find_if(std::begin(__s), std::end(__s),
                                [](const double d) { return d != 0.; })
                   == std::end(__s))),
              s(std::move(__s)), h(std::move(__h)), g(std::move(gen))
        {
        }
        bool
        operator==(const generalmut &rhs) const
        {
            return this->pos == rhs.pos && this->neutral == rhs.neutral
                   && this->s == rhs.s && this->h == rhs.h;
        }
    };

    /*!
      \brief Mutation type allowing arbitray number of "s,h" pairs representing
      effect sizes,dominance.
      \ingroup sugar

      This type represents a mutation with a number of 's' and 'h' parameters
      that are determined at run time.

      This type complements KTfwd::generalmut.  The use cases will be similar.
      This type will show better performance
      for large numbers of categories (e.g., large sizes of member variables s
      and h).

      This type differs from KTfwd::generalmut in that the sizes of the s and h
      vectors do not have to be the same.
     */
    struct generalmut_vec : public mutation_base
    {
        using array_t = std::vector<double>;
        //! Selection coefficients and/or effect sizes
        array_t s;
        //! Dominances associated w/values in 's'.
        array_t h;
        //! Generation when mutation arose
        uint_t g;
        //! Constructor
        generalmut_vec(array_t &&__s, array_t &&__h, double pos, uint_t gen)
            : KTfwd::mutation_base(
                  std::move(pos),
                  // Mutation is neutral i.f.f. all values in __s == 0.
                  (std::find_if(std::begin(__s), std::end(__s),
                                [](const double d) { return d != 0.; })
                   == std::end(__s))),
              s(std::move(__s)), h(std::move(__h)), g(std::move(gen))
        {
        }
        bool
        operator==(const generalmut_vec &rhs) const
        {
            return this->pos == rhs.pos && this->neutral == rhs.neutral
                   && this->s == rhs.s && this->h == rhs.h;
        }
    };
}
#endif
