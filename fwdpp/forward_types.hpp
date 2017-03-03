/*!
  \file forward_types.hpp
  \defgroup basicTypes Mutation and gamete data types
  These are the basic types common to all simulations: mutations and gametes.
  See @ref md_md_policies for more details.
*/
#ifndef _FORWARD_TYPES_HPP_
#define _FORWARD_TYPES_HPP_

#include <limits>
#include <vector>
#include <cmath>
#include <cstdint>
#include <type_traits>
#include <fwdpp/tags/tags.hpp>

namespace KTfwd
{
    //! The unsigned integer type is 32 bits
    using uint_t = std::uint32_t;

    /*! \brief Base class for mutations
      At minimum, a mutation must contain a position and a count in the
      population.
      You can derive from this class, for instance to add selection
      coefficients,
      counts in different sexes, etc.
      \ingroup basicTypes
      \note See @ref TutMut in @ref md_md_policies for more detail on how to
      extend this type
    */
    struct mutation_base
    {
        /// Mutation position
        double pos;
        /*!
          16 bits of extra data to be associated w/this type.
          Do with it what you will. Fits into padded space in this struct,
          and doesn't affect sizeof(mutation).
        */
        std::uint16_t xtra;
        /// Is the mutation neutral or not?
        bool neutral;
        mutation_base(const double &position, const bool &isneutral = true,
                      const std::uint16_t x = 0) noexcept : pos(position),
                                                            xtra(x),
                                                            neutral(isneutral)
        {
        }
        virtual ~mutation_base() noexcept {}
        mutation_base(mutation_base &) = default;
        mutation_base(mutation_base const &) = default;
        mutation_base(mutation_base &&) = default;
        mutation_base &operator=(mutation_base &) = default;
        mutation_base &operator=(mutation_base const &) = default;
        mutation_base &operator=(mutation_base &&) = default;
    };

    struct mutation : public mutation_base
    /*!
      \brief The simplest mutation type, adding just a selection coefficient
      and dominance to the interface
      \ingroup basicTypes
    */
    {
        /// selection coefficient
        double s;
        /// dominance coefficient
        double h;
        mutation(const double &position, const double &sel_coeff,
                 const double &dominance = 0.5) noexcept
            : mutation_base(position, (sel_coeff == 0)),
              s(sel_coeff),
              h(dominance)
        {
        }
        bool
        operator==(const mutation &rhs) const
        {
            return (std::fabs(this->pos - rhs.pos)
                        <= std::numeric_limits<double>::epsilon()
                    && this->s == rhs.s);
        }
    };

    /*! \brief Base class for gametes.

      A gamete contains one container of keys to neutral mutations, and another
          to selected mutations.  It also keeps track of its count (number of
      occurrences).

      The template parameter types are:
      tag_type = A type that can be used as a "dispatch tag".  Currently, these
      are not used elsewhere in the library, but they may
      be in the future, or this may disappear in future library releases.

      \note The typical use of this class is simply to define your mutation
      type (see @ref md_md_policies)
      and then use a typedef to define your gamete type in the simulations:
      \code
      using gamete_t = KTfwd::gamete_base<mutation_type>
      \endcode
      See @ref md_md_policies for examples of this.
      \ingroup basicTypes
    */
    template <typename TAG = tags::standard_gamete> struct gamete_base
    {
        //! Count in population
        uint_t n;
        //! Dispatch tag type
        using gamete_tag = TAG;
        using index_t = std::uint32_t;
        using mutation_container = std::vector<index_t>;
        //! Container of mutations not affecting trait value/fitness ("neutral
        //! mutations")
        mutation_container mutations;
        //! Container of mutations affecting trait value/fitness ("selected"
        //! mutations")
        mutation_container smutations;

        /*! @brief Constructor
          \param icount The number of occurrences of this gamete in the
          population
        */
        gamete_base(const uint_t &icount) noexcept
            : n(icount),
              mutations(mutation_container()),
              smutations(mutation_container())
        {
        }

        /*! @brief "Perfect-forwarding" constructor
          \param icount The number of occurrences of this gamete in the
          population
          \param n A container of mutations not affecting trait value/fitness
          \param s A container of mutations affecting trait value/fitness
        */
        template <typename T>
        gamete_base(const uint_t &icount, T &&n, T &&s) noexcept
            : n(icount),
              mutations(std::forward<T>(n)),
              smutations(std::forward<T>(s))
        {
        }
        //! Destructor is virtual, so you may inherit from this type
        virtual ~gamete_base() noexcept {}
        //! Copy constructor
        gamete_base(gamete_base &) = default;
        //! Copy constructor
        gamete_base(gamete_base const &) = default;
        //! Move constructor
        gamete_base(gamete_base &&) = default;

        //! Assignment operator
        gamete_base &operator=(gamete_base &) = default;
        //! Assignment operator
        gamete_base &operator=(gamete_base const &) = default;
        //! Move assignment operator
        gamete_base &operator=(gamete_base &&) = default;
        /*! \brief Equality operation
        */
        inline bool
        operator==(const gamete_base<TAG> &rhs) const
        {
            return (this->mutations == rhs.mutations
                    && this->smutations == rhs.smutations);
        }
// The following type traits are not implemented in GCC 4.
#if (defined __GNUG__&& __GNUC__ >= 5) || !defined __GNUG__ || defined __clang__ 
        static_assert(std::is_trivially_constructible<gamete_tag>::value,
                      "Typename gamete_tag must refer to a "
                      "trivially-constrictible type");
        static_assert(std::is_nothrow_default_constructible<gamete_tag>::value,
                      "Typename gamete_tag must refer to a type that is "
                      "default- and nothrow-constructible");
#endif
    };

    /// Default gamete type
    using gamete = gamete_base<tags::standard_gamete>;
}
#endif /* _FORWARD_TYPES_HPP_ */
