#ifndef __FWDPP_SUGAR_MUTATION_POPGENMUT_HPP__
#define __FWDPP_SUGAR_MUTATION_POPGENMUT_HPP__

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/tags/tags.hpp>
#include <fwdpp/io/scalar_serialization.hpp>
#include <fwdpp/io/mutation.hpp>
#include <fwdpp/simfunctions/recycling.hpp>
#include <limits>
#include <tuple>

namespace fwdpp
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
        //! Alias for tuple type that can be used for object construction
        using constructor_tuple
            = std::tuple<double, double, double, unsigned, std::uint16_t>;

        /*!
          \brief Constructor
          \param __pos Mutation position
          \param __s Selection coefficient
          \param __h Dominance coefficient
          \param __g Generation when mutation arose
          \param x Value to assign to mutation_base::xtra
        */
        popgenmut(const double &__pos, const double &__s, const double &__h,
                  const unsigned &__g, const std::uint16_t x = 0) noexcept
            : mutation_base(__pos, (__s == 0.) ? true : false, x), g(__g),
              s(__s), h(__h)
        {
        }

        ///
        /// Construct from a tuple.
        ///
        /// \param t Elements must be pos, s, h, g, x
        ///
        /// \version
        /// Added in fwdpp 0.5.7
        popgenmut(constructor_tuple t) noexcept
            : mutation_base(std::get<0>(t),
                            (std::get<1>(t) == 0.) ? true : false,
                            std::get<4>(t)),
              g(std::get<3>(t)), s(std::get<1>(t)), h(std::get<2>(t))
        {
        }

        bool
        operator==(const popgenmut &rhs) const
        {
            return std::tie(this->g, this->s, this->h)
                       == std::tie(rhs.g, rhs.s, rhs.h)
                   && is_equal(rhs);
        }
    };

    template < typename mcont_t, typename lookup_table_t,
              typename position_function, typename effect_size_function,
              typename dominance_function>
    std::size_t
    infsites_popgenmut(flagged_mutation_queue &recycling_bin, mcont_t &mutations,
                       const gsl_rng *r, lookup_table_t &lookup,
                       const uint_t &generation, const double pselected,
                       const position_function &posmaker,
                       const effect_size_function &esize_maker,
                       const dominance_function &hmaker,
                       const decltype(popgenmut::xtra) x = 0)
    /*!
	 * Mutation function to add a fwdpp::popgenmut to a population.
	 *
	 * In order to use this function, it must be bound to a callable
	 * that is a valid mutation function.  See examples for details.
	 *
	 * \param recycling_bin Recycling queue for mutations (fwdpp::flagged_mutation_queue).
	 * \param mutations Container of mutations
	 * \param r A random-number generator
	 * \param lookup Lookup table for mutation positions
	 * \param generation The generation that is being mutated
	 * \param pselected  The probability that a new mutation affects fitness
	 * \param posmaker A function generating a mutation position.  Must be convertible to std::function<double()>.
	 * \param esize_maker A function to generate an effect size, given that a mutation affects fitness. Must be convertible to std::function<double()>.
	 * \param hmaker A function to generate a dominance value, given that a mutation affects fitness. Must be convertible to std::function<double()>.
     * \param x Value to pass as popgenmut::xtra
	 *
	 * \note "Neutral" mutations get assigned a dominance of zero.  The xtra field is not written to.
	 *
	 * \version 0.6.0
	 * Added to library
	 */
    {
        auto pos = posmaker();
        while (lookup.find(pos) != lookup.end())
            {
                pos = posmaker();
            }
        bool selected = (gsl_rng_uniform(r) < pselected);
        auto idx = recycle_mutation_helper(
            recycling_bin, mutations, pos, (selected) ? esize_maker() : 0.,
            (selected) ? hmaker() : 0., generation, x);
        lookup.emplace(pos, idx);
        return idx;
    }

    namespace io
    {
        template <> struct serialize_mutation<popgenmut>
        /// Specialization for fwdpp::popgenmut
        {
            io::scalar_writer writer;
            serialize_mutation<popgenmut>() : writer{} {}
            template <typename streamtype>
            inline void
            operator()(streamtype &buffer, const popgenmut &m) const
            {
                writer(buffer, &m.g);
                writer(buffer, &m.pos);
                writer(buffer, &m.s);
                writer(buffer, &m.h);
                writer(buffer, &m.xtra);
            }
        };

        template <> struct deserialize_mutation<popgenmut>
        /// Specialization for fwdpp::popgenmut
        {
            io::scalar_reader reader;
            deserialize_mutation<popgenmut>() : reader{} {}
            template <typename streamtype>
            inline popgenmut
            operator()(streamtype &buffer) const
            {
                uint_t g;
                double pos, s, h;
                decltype(popgenmut::xtra) xtra;
                io::scalar_reader reader;
                reader(buffer, &g);
                reader(buffer, &pos);
                reader(buffer, &s);
                reader(buffer, &h);
                reader(buffer, &xtra);

                return popgenmut(pos, s, h, g, xtra);
            }
        };
    } // namespace io
} // namespace fwdpp
#endif
