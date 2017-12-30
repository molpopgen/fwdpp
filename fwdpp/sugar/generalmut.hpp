#ifndef FWDPP_SUGAR_GENERALMUT_HPP
#define FWDPP_SUGAR_GENERALMUT_HPP

#include <tuple>
#include <array>
#include <algorithm>
#include <memory>
#include <limits>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/io/scalar_serialization.hpp>
#include <fwdpp/io/mutation.hpp>
#include <fwdpp/tags/tags.hpp>

namespace fwdpp
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
      at run-time.  For that use case, see fwdpp::generalmut_vec, which uses
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
        //! Tuple type useable for object construction
        using constructor_tuple
            = std::tuple<array_t, array_t, double, uint_t, std::uint16_t>;

        //! Constructor
        generalmut(array_t __s, array_t __h, double pos, uint_t gen,
                   std::uint16_t label = 0)
            : mutation_base(
                  pos,
                  // Mutation is neutral i.f.f. all values in __s == 0.
                  (std::find_if(std::begin(__s), std::end(__s),
                                [](const double d) { return d != 0.; })
                   == std::end(__s)),
                  label),
              s(std::move(__s)), h(std::move(__h)), g(std::move(gen))
        {
        }

        ///
        /// Constructor from a tuple
        ///
        /// \param t A tuple (s,g,pos,origin time, xtra)
        ///
        /// \version
        /// Added in fwdpp 0.5.7
        generalmut(constructor_tuple t)
            : mutation_base(
                  std::get<2>(t),
                  (std::find_if(std::begin(std::get<0>(t)),
                                std::end(std::get<0>(t)),
                                [](const double d) { return d != 0.; })
                   == std::end(std::get<0>(t))),
                  std::get<4>(t)),
              s(std::move(std::get<0>(t))), h(std::move(std::get<1>(t))),
              g(std::get<3>(t))
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

      This type complements fwdpp::generalmut.  The use cases will be similar.
      This type will show better performance
      for large numbers of categories (e.g., large sizes of member variables s
      and h).

      This type differs from fwdpp::generalmut in that the sizes of the s and h
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
        //! Tuple type useable for object construction
        using constructor_tuple
            = std::tuple<array_t, array_t, double, uint_t, std::uint16_t>;
        //! Constructor
        generalmut_vec(array_t &&__s, array_t &&__h, double pos, uint_t gen,
                       const std::uint16_t x = 0)
            : fwdpp::mutation_base(
                  std::move(pos),
                  // Mutation is neutral i.f.f. all values in __s == 0.
                  (std::find_if(std::begin(__s), std::end(__s),
                                [](const double d) { return d != 0.; })
                   == std::end(__s)),
                  x),
              s(std::move(__s)), h(std::move(__h)), g(std::move(gen))
        {
        }

        ///
        /// Constructor from a tuple
        ///
        /// \param t A tuple (s,g,pos,origin time, xtra)
        ///
        /// \version
        /// Added in fwdpp 0.5.7
        generalmut_vec(constructor_tuple t)
            : mutation_base(
                  std::get<2>(t),
                  (std::find_if(std::begin(std::get<0>(t)),
                                std::end(std::get<0>(t)),
                                [](const double d) { return d != 0.; })
                   == std::end(std::get<0>(t))),
                  std::get<4>(t)),
              s(std::move(std::get<0>(t))), h(std::move(std::get<1>(t))),
              g(std::get<3>(t))
        {
        }

        bool
        operator==(const generalmut_vec &rhs) const
        {
            return this->pos == rhs.pos && this->neutral == rhs.neutral
                   && this->s == rhs.s && this->h == rhs.h;
        }
    };

    namespace io
    {
        template <> struct serialize_mutation<generalmut_vec>
        {
            io::scalar_writer writer;
            serialize_mutation<generalmut_vec>() : writer{} {}
            template <typename streamtype>
            inline void
            operator()(const generalmut_vec &m, streamtype &buffer) const
            {
                writer(buffer, &m.pos);
                writer(buffer, &m.g);
                writer(buffer, &m.xtra);
                // Write mutation types
                using array_t_size_t =
                    typename generalmut_vec::array_t::size_type;
                array_t_size_t ns = m.s.size(), nh = m.h.size();
                writer(buffer, &ns);
                writer(buffer, &nh);
                writer(buffer, m.s.data(), ns);
                writer(buffer, m.h.data(), nh);
            }
        };

        template <> struct deserialize_mutation<generalmut_vec>
        {
            io::scalar_reader reader;
            deserialize_mutation<generalmut_vec>() : reader{} {}
            template <typename streamtype>
            inline generalmut_vec
            operator()(streamtype &buffer) const
            {
                double pos;
                decltype(generalmut_vec::g) g;
                decltype(generalmut_vec::xtra) xtra;
                io::scalar_reader reader;
                reader(buffer, &pos);
                reader(buffer, &g);
                reader(buffer, &xtra);
                typename generalmut_vec::array_t::size_type ns, nh;
                reader(buffer, &ns);
                reader(buffer, &nh);
                typename generalmut_vec::array_t s(ns), h(nh);
                // Write mutation types
                if (ns)
                    {
                        reader(buffer, s.data(), ns);
                    }
                if (nh)
                    {
                        reader(buffer, h.data(), nh);
                    }
                return generalmut_vec(std::move(s), std::move(h), pos, g,
                                      xtra);
            }
        };
    }
}
#define SPECIALIZE_SERIALIZE_MUTATION_GENERALMUT_BODY(N)                      \
    template <typename streamtype>                                            \
    inline void operator()(const generalmut<N> &m, streamtype &buffer) const  \
    {                                                                         \
        io::scalar_writer writer;                                             \
        writer(buffer, &m.g);                                                 \
        writer(buffer, &m.pos);                                               \
        writer(buffer, &m.xtra);                                              \
        writer(buffer, &m.s[0], N);                                           \
        writer(buffer, &m.h[0], N);                                           \
    }

#define SPECIALIZE_DESERIALIZE_MUTATION_GENERALMUT_BODY(N)                    \
    template <typename streamtype>                                            \
    inline generalmut<N> operator()(streamtype &buffer) const                 \
    {                                                                         \
        uint_t g;                                                             \
        double pos;                                                           \
        decltype(generalmut<N>::xtra) xtra;                                   \
        using value_t = generalmut<N>::array_t::value_type;                   \
        io::scalar_reader reader;                                             \
        std::array<value_t, std::tuple_size<generalmut<N>::array_t>::value>   \
            s, h;                                                             \
        reader(buffer, &g);                                                   \
        reader(buffer, &pos);                                                 \
        reader(buffer, &xtra);                                                \
        reader(buffer, &s[0],                                                 \
               std::tuple_size<generalmut<N>::array_t>::value);               \
        reader(buffer, &h[0],                                                 \
               std::tuple_size<generalmut<N>::array_t>::value);    \
        return generalmut<N>(s, h, pos, g, xtra);                             \
    }

#endif
