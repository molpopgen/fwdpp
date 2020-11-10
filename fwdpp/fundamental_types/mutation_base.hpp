#ifndef FWDPP_FUNDAMENTAL_TYPES_MUTATION_BASE_HPP__
#define FWDPP_FUNDAMENTAL_TYPES_MUTATION_BASE_HPP__

#include <cstdint>

namespace fwdpp
{
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
                      const std::uint16_t x = 0) noexcept
            : pos(position), xtra(x), neutral(isneutral)
        {
        }
        virtual ~mutation_base() noexcept
        {
        }
        mutation_base(mutation_base const &) = default;
        mutation_base(mutation_base &&) = default;
        mutation_base &operator=(mutation_base const &) = default;
        mutation_base &operator=(mutation_base &&) = default;

        inline bool
        is_equal(const mutation_base &rhs) const
        {
            return this->pos == rhs.pos && this->xtra == rhs.xtra
                   && this->neutral == rhs.neutral;
        }
    };
}

#endif
