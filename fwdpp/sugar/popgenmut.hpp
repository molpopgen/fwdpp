#ifndef __FWDPP_SUGAR_MUTATION_POPGENMUT_HPP__
#define __FWDPP_SUGAR_MUTATION_POPGENMUT_HPP__

#include <gsl/gsl_rng.h>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/tags/tags.hpp>
#include <fwdpp/io/scalar_serialization.hpp>
#include <fwdpp/io/mutation.hpp>
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
          \param __x Value to assign to mutation_base::xtra
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

        template <typename mcont_t, typename lookup_table_t, typename posmaker,
                  typename smaker, typename hmaker>
        inline static std::function<std::size_t(
            fwdpp::traits::recycling_bin_t<mcont_t> &, mcont_t &)>
        infsites(const gsl_rng *r, double psel, const mcont_t &,
                 lookup_table_t &lookup, uint_t &generation, posmaker p,
                 smaker s, hmaker h)
        {
            auto mutation_model = [&r, &generation, &lookup, psel, p, s, h](
                fwdpp::traits::recycling_bin_t<mcont_t> &recycling_bin,
                mcont_t &mutations) {
                auto pos = p();
                while (lookup.find(pos) != lookup.end())
                    {
                        pos = p();
                    }
                lookup.insert(pos);
                bool selected = gsl_rng_uniform(r) < psel;
                return fwdpp_internal::recycle_mutation_helper(
                    recycling_bin, mutations, pos, (selected) ? s() : 0.,
                    (selected) ? h() : 0., generation);
            };
            return mutation_model;
        }
    };

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
    }
}
#endif
