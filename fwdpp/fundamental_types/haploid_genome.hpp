#ifndef FWDPP_FUNDAMENTAL_TYPES_HAPLOID_GENOME_HPP__
#define FWDPP_FUNDAMENTAL_TYPES_HAPLOID_GENOME_HPP__

#include <tuple>
#include <vector>
#include <type_traits>
#include <fwdpp/tags/tags.hpp>
#include "typedefs.hpp"

namespace fwdpp
{
    /*! \brief Base class for genomes.

      A haploid_genome contains one container of keys to neutral mutations, and another
          to selected mutations.  It also keeps track of its count (number of
      occurrences).

      The template parameter types are:
      tag_type = A type that can be used as a "dispatch tag".  Currently, these
      are not used elsewhere in the library, but they may
      be in the future, or this may disappear in future library releases.

      \ingroup basicTypes
    */
    template <typename TAG = tags::standard_haploid_genome> struct haploid_genome_base
    {
        //! Count in population
        uint_t n;
        //! Dispatch tag type
        using haploid_genome_tag = TAG;
        /// Container type for mutation indexes
        using mutation_container = std::vector<uint_t>;
        /// The integer type used to store mutation indexes
        using index_t = typename mutation_container::value_type;
        //! Container of mutations not affecting trait value/fitness
        mutation_container mutations;
        //! Container of mutations affecting trait value/fitness
        mutation_container smutations;

        //! Tuple type usable for construction
        using constructor_tuple
            = std::tuple<uint_t, mutation_container, mutation_container>;

        /*! @brief Constructor
          \param icount The number of occurrences of this haploid_genome in the
          population
        */
        explicit haploid_genome_base(const uint_t &icount) noexcept
            : n(icount), mutations(mutation_container()),
              smutations(mutation_container())
        {
        }

        /*! @brief "Perfect-forwarding" constructor
          \param icount The number of occurrences of this haploid_genome in the
          population
          \param n A container of mutations not affecting trait value/fitness
          \param s A container of mutations affecting trait value/fitness
        */
        template <typename T>
        haploid_genome_base(const uint_t &icount, T &&n, T &&s) noexcept
            : n(icount), mutations(std::forward<T>(n)), smutations(std::forward<T>(s))
        {
        }

        //! Construct from n, mutations, smutations tuple
        [[deprecated]] haploid_genome_base(constructor_tuple t)
            : n(std::get<0>(t)), mutations(std::move(std::get<1>(t))),
              smutations(std::move(std::get<2>(t)))
        {
        }

        //! Destructor is virtual, so you may inherit from this type
        virtual ~haploid_genome_base() noexcept
        {
        }

        //! Copy constructor
        haploid_genome_base(haploid_genome_base const &) = default;
        //! Move constructor
        haploid_genome_base(haploid_genome_base &&) = default;

        //! Assignment operator
        haploid_genome_base &operator=(haploid_genome_base const &) = default;
        //! Move assignment operator
        haploid_genome_base &operator=(haploid_genome_base &&) = default;
        /*! \brief Equality operation
        */
        inline bool
        operator==(const haploid_genome_base<TAG> &rhs) const
        {
            return (this->mutations == rhs.mutations
                    && this->smutations == rhs.smutations);
        }
// The following type traits are not implemented in GCC 4.
#if (defined __GNUG__ && __GNUC__ >= 5) || !defined __GNUG__ || defined __clang__
        static_assert(std::is_trivially_constructible<haploid_genome_tag>::value,
                      "Typename haploid_genome_tag must refer to a "
                      "trivially-constrictible type");
        static_assert(std::is_nothrow_default_constructible<haploid_genome_tag>::value,
                      "Typename haploid_genome_tag must refer to a type that is "
                      "default- and nothrow-constructible");
#endif
    };

    /// \typedef haploid_genome
    /// Default haploid_genome type
    using haploid_genome = haploid_genome_base<tags::standard_haploid_genome>;
}

#endif
