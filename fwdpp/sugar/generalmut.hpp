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
        using array_t = std::array<std::tuple<double, double>, N>;
        //! Effect sizes and dominance
        array_t sh;
        //! Generation when mutation arose
        uint_t g;
        //! Tuple type useable for object construction
        using constructor_tuple
            = std::tuple<array_t, double, uint_t, std::uint16_t>;

        //! Constructor
        generalmut(array_t __sh, double pos, uint_t gen,
                   std::uint16_t label = 0)
            : mutation_base(
                  pos,
                  // Mutation is neutral i.f.f. all values in __s == 0.
                  (std::find_if(std::begin(__sh), std::end(__sh),
                                [](const typename array_t::value_type &t) {
                                    return std::get<0>(t) != 0.;
                                })
                   == std::end(__sh)),
                  label),
              sh(std::move(__sh)), g{ gen }
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
                  std::get<1>(t),
                  (std::find_if(std::begin(std::get<0>(t)),
                                std::end(std::get<0>(t)),
                                [](const std::tuple<double, double> &t) {
                                    return std::get<0>(t) != 0.;
                                })
                   == std::end(std::get<0>(t))),
                  std::get<3>(t)),
              sh(std::move(std::get<0>(t))), g(std::get<2>(t))
        {
        }

        bool
        operator==(const generalmut &rhs) const
        {
            return this->g == rhs.g && this->sh == rhs.sh && is_equal(rhs);
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
        using array_t = std::vector<std::tuple<double, double>>;
        //! Effect size and dominance tuples:
        array_t sh;
        //! Generation when mutation arose
        uint_t g;
        //! Tuple type useable for object construction
        using constructor_tuple
            = std::tuple<array_t, double, uint_t, std::uint16_t>;
        //! Constructor
        generalmut_vec(array_t __sh, double pos, uint_t gen,
                       const std::uint16_t x = 0)
            : fwdpp::mutation_base(
                  std::move(pos),
                  // Mutation is neutral i.f.f. all values in __s == 0.
                  (std::find_if(std::begin(__sh), std::end(__sh),
                                [](const array_t::value_type t) {
                                    return std::get<0>(t) != 0.;
                                })
                   == std::end(__sh)),
                  x),
              sh(std::move(__sh)), g(std::move(gen))
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
            : mutation_base(std::get<1>(t),
                            (std::find_if(std::begin(std::get<0>(t)),
                                          std::end(std::get<0>(t)),
                                          [](const array_t::value_type t) {
                                              return std::get<0>(t) != 0.;
                                          })
                             == std::end(std::get<0>(t))),
                            std::get<3>(t)),
              sh(std::move(std::get<0>(t))), g(std::get<2>(t))
        {
        }

        bool
        operator==(const generalmut_vec &rhs) const
        {
            return this->g == rhs.g && this->sh == rhs.sh && is_equal(rhs);
        }
    };

    namespace io
    {
        template <> struct serialize_mutation<generalmut_vec>
        /// Specialization for fwdpp::generalmut_vec
        {
            io::scalar_writer writer;
            serialize_mutation<generalmut_vec>() : writer{} {}
            template <typename streamtype>
            inline void
            operator()(streamtype &buffer, const generalmut_vec &m) const
            {
                writer(buffer, &m.pos);
                writer(buffer, &m.g);
                writer(buffer, &m.xtra);
                // Write mutation types
                using array_t_size_t =
                    typename generalmut_vec::array_t::size_type;
                array_t_size_t ns = m.sh.size();
                writer(buffer, &ns);
                for (auto &&sh : m.sh)
                    {
                        writer(buffer, &std::get<0>(sh));
                        writer(buffer, &std::get<1>(sh));
                    }
            }
        };

        template <> struct deserialize_mutation<generalmut_vec>
        /// Specialization for fwdpp::generalmut_vec
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
                typename generalmut_vec::array_t::size_type ns;
                reader(buffer, &ns);
                typename generalmut_vec::array_t sh(ns);
                // Write mutation types
                if (ns)
                    {
                        for (auto &i : sh)
                            {
                                reader(buffer, &std::get<0>(i));
                                reader(buffer, &std::get<1>(i));
                            }
                    }
                return generalmut_vec(std::move(sh), pos, g, xtra);
            }
        };
    }
}
#define SPECIALIZE_SERIALIZE_MUTATION_GENERALMUT_BODY(N)                      \
    template <typename streamtype>                                            \
    inline void operator()(streamtype &buffer, const generalmut<N> &m) const  \
    {                                                                         \
        io::scalar_writer writer;                                             \
        writer(buffer, &m.g);                                                 \
        writer(buffer, &m.pos);                                               \
        writer(buffer, &m.xtra);                                              \
        for (auto &&sh : m.sh)                                                \
            {                                                                 \
                writer(buffer, &std::get<0>(sh));                             \
                writer(buffer, &std::get<1>(sh));                             \
            }                                                                 \
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
            sh;                                                               \
        reader(buffer, &g);                                                   \
        reader(buffer, &pos);                                                 \
        reader(buffer, &xtra);                                                \
        for (auto &sh_i : sh)                                                 \
            {                                                                 \
                reader(buffer, &std::get<0>(sh_i));                           \
                                                                              \
                reader(buffer, &std::get<1>(sh_i));                           \
            }                                                                 \
        return generalmut<N>(sh, pos, g, xtra);                               \
    }

#endif
