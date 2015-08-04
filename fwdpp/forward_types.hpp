/*!
  \file forward_types.hpp
  \defgroup basicTypes Mutation and gamete data types
  These are the basic types common to all simulations: mutations and gametes.  See @ref md_md_policies for more details.
*/
#ifndef _FORWARD_TYPES_HPP_
#define _FORWARD_TYPES_HPP_

#include <limits>
#include <vector>
#include <list>
#include <cmath>
#include <type_traits>
#include <fwdpp/tags/gamete_tags.hpp>

namespace KTfwd
{
  /*! \brief Base class for mutations
    At minimum, a mutation must contain a position and a count in the population.	
    You can derive from this class, for instance to add selection coefficients,
    counts in different sexes, etc.
    \ingroup basicTypes
    \note See @ref TutMut in @ref md_md_policies for more detail on how to extend this type
  */
  struct mutation_base
  {
    /// Mutation position
    mutable double pos;
    /// Count of mutation in the population
    unsigned n;
    /// Is the mutation neutral or not?
    bool neutral;
    /// Used internally (don't worry about it for now...)
    bool checked;
    mutation_base(const double & position, const unsigned & count, const bool & isneutral = true) noexcept
      : pos(position),n(count),neutral(isneutral),checked(false)
    {	
    }
    virtual ~mutation_base() noexcept {}
    mutation_base( mutation_base & ) = default;
    mutation_base( mutation_base const & ) = default;
    mutation_base( mutation_base && ) = default;
    mutation_base & operator=(mutation_base &) = default;
    mutation_base & operator=(mutation_base const &) = default;
    mutation_base & operator=(mutation_base &&) = default;
  };

  struct mutation : public mutation_base
  /*!
    \brief The simplest mutation type, adding just a selection coefficient and dominance to the interface
    \ingroup basicTypes
  */
  {
    /// selection coefficient
    mutable double s;
    /// dominance coefficient
    mutable double h;
    mutation( const double & position, const double & sel_coeff,const unsigned & count,
	      const double & dominance = 0.5) 
      : mutation_base(position,count,(sel_coeff==0)),s(sel_coeff),h(dominance)
    {
    }
    bool operator==(const mutation & rhs) const
    {
      return( std::fabs(this->pos-rhs.pos) <= std::numeric_limits<double>::epsilon() &&
	      this->s == rhs.s );
    }
  };

  /*! \brief Base class for gametes.

    A gamete is a container of pointers (iterators) to mutations + a count in the population.

    The template parameter types are:
    mut_type = the mutation type to be used.  Must be a model of KTfwd::mutation_base
    list_type = the (doubly-linked) list that mutations are stored in.  This is mainly used for defining types for this class
    tag_type = A type that can be used as a "dispatch tag".  Currently, these are not used elsewhere in the library, but they may
    be in the future, or this may disappear in future library releases.  The current default (KTfwd::tags::standard_gamete) maintains
    backwards compatibility with previous library versions and does not affect compilation of existing programs based on the library.

    \note The typical use of this class is simply to define your mutation type (see @ref md_md_policies)
    and then use a typedef to define your gamete type in the simulations:
    \code
    using gamete_t = KTfwd::gamete_base<mutation_type>
    \endcode
    See @ref md_md_policies for examples of this.
    \ingroup basicTypes
  */
  template<typename mut_type,
	   typename list_type = std::list<mut_type>,
	   typename tag_type = tags::standard_gamete>
  struct gamete_base
  {
    static_assert( std::is_base_of<mutation_base,mut_type>::value,
		   "mut_type must be derived from KTfwd::mutation_base" );
    //! Count in population
    unsigned n;
    using mutation_type = mut_type;
    using mutation_list_type = list_type;
    using mutation_list_type_iterator = typename list_type::iterator;
    using mutation_container = std::vector< mutation_list_type_iterator >;
    using mcont_iterator = typename mutation_container::iterator;
    using mcont_const_iterator = typename mutation_container::const_iterator;
    //! Dispatch tag type
    using gamete_tag = tag_type;
    //! Container of mutations not affecting trait value/fitness ("neutral mutations")
    mutation_container mutations;
    //! Container of mutations affecting trait value/fitness ("neutral mutations")
    mutation_container smutations;

    /*! @brief Constructor
      \param icount The number of occurrences of this gamete in the population
    */
    gamete_base(const unsigned & icount) noexcept : n(icount),mutations( mutation_container() ),smutations(mutation_container())
    {
    }

    /*! @brief Constructor
      \param icount The number of occurrences of this gamete in the population
      \param n A container of mutations not affecting trait value/fitness
      \param s A container of mutations affecting trait value/fitness
    */
    gamete_base(const unsigned & icount, const mutation_container & n,
		const mutation_container & s) noexcept : n(icount),mutations(n),smutations(s)
    {
    }
    //! Destructor is virtual, so you may inherit from this type
    virtual ~gamete_base() noexcept {}
    //! Copy constructor
    gamete_base( gamete_base & ) = default;
    //! Copy constructor
    gamete_base( gamete_base const & ) = default;
    //! Move constructor
    gamete_base( gamete_base && ) = default;

    //! Assignment operator
    gamete_base & operator=(gamete_base &) = default;
    //! Assignment operator
    gamete_base & operator=(gamete_base const &) = default;
    //! Move assignment operator
    gamete_base & operator=(gamete_base &&) = default;
    /*! \brief Equality operation
      \note Given that mutations and smutations contains ITERATORS to actual mutations,
      operator== does not need to be defined for the corresponding mutation type
    */
    inline bool operator==(const gamete_base<mut_type,list_type> & rhs) const
    {
      return(this->mutations == rhs.mutations && this->smutations == rhs.smutations);
    }
  };

  /*! The simplest gamete adds nothing to the interface of the base class.
    \ingroup basicTypes
  */
  using gamete = gamete_base<mutation>;

}
#endif /* _FORWARD_TYPES_HPP_ */
