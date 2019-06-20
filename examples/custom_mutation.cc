/// \include custom_mutation.cc
/// An example of a custom mutation type,
/// including serialization and
/// a mutation generating function.
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <array>
#include <algorithm>
#include <limits>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/io/scalar_serialization.hpp>
#include <fwdpp/io/mutation.hpp>
#include <fwdpp/simfunctions/recycling.hpp>
#include <fwdpp/type_traits.hpp>

struct TwoDmutation : public fwdpp::mutation_base
{
    using array_t = std::array<double, 2>;
    //! Effect sizes and dominance
    array_t s, h;
    //! Generation when mutation arose
    fwdpp::uint_t g;

    //! Constructor
    TwoDmutation(array_t s_, array_t h_, double pos, fwdpp::uint_t gen,
                 std::uint16_t label = 0)
        : mutation_base(pos,
                        // Mutation is neutral i.f.f. all values in s_ == 0.
                        !std::any_of(std::begin(s_), std::end(s_),
                                     [](const double t) { return t != 0.; }),
                        label),
          s(std::move(s_)), h(std::move(h_)), g{ gen }
    {
    }

    bool
    operator==(const TwoDmutation &rhs) const
    {
        return this->g == rhs.g && this->s == rhs.s && this->h == rhs.h
               && is_equal(rhs);
    }
    bool
    operator!=(const TwoDmutation &rhs) const
    {
        return !(*this == rhs);
    }
};

namespace fwdpp
{
    namespace io
    {
        template <> struct serialize_mutation<TwoDmutation>
        /// Specialization for fwdpp::TwoDmutation
        {
            io::scalar_writer writer;
            serialize_mutation<TwoDmutation>() : writer{} {}
            template <typename streamtype>
            inline void
            operator()(streamtype &buffer, const TwoDmutation &m) const
            {
                writer(buffer, &m.pos);
                writer(buffer, &m.g);
                writer(buffer, &m.xtra);
                // Write mutation data
                writer(buffer, m.s.data(), 2);
                writer(buffer, m.h.data(), 2);
            }
        };

        template <> struct deserialize_mutation<TwoDmutation>
        /// Specialization for fwdpp::TwoDmutation
        {
            io::scalar_reader reader;
            deserialize_mutation<TwoDmutation>() : reader{} {}
            template <typename streamtype>
            inline TwoDmutation
            operator()(streamtype &buffer) const
            {
                double pos;
                decltype(TwoDmutation::g) g;
                decltype(TwoDmutation::xtra) xtra;
                io::scalar_reader reader;
                reader(buffer, &pos);
                reader(buffer, &g);
                reader(buffer, &xtra);
                TwoDmutation::array_t s, h;
                reader(buffer, s.data(), 2);
                reader(buffer, h.data(), 2);
                return TwoDmutation(std::move(s), std::move(h), pos, g, xtra);
            }
        };
    } // namespace io
} // namespace fwdpp

template <typename mcont_t, typename lookup_table_t,
          typename position_function>
std::size_t
infsites_TwoDmutation_additive(
    fwdpp::flagged_mutation_queue &recycling_bin, mcont_t &mutations,
    const gsl_rng *r, lookup_table_t &lookup, const fwdpp::uint_t &generation,
    const double pselected, const position_function &posmaker,
    const double sigma_x, const double sigma_y, const double rho,
    const decltype(TwoDmutation::xtra) x = 0)
{
    TwoDmutation::array_t s{ 0., 0. }, h{ 1., 1. };

    double pos = posmaker();
    while (lookup.find(pos) != lookup.end())
        {
            pos = posmaker();
        }

    if (gsl_rng_uniform(r) < pselected)
        {
            gsl_ran_bivariate_gaussian(r, sigma_x, sigma_y, rho, &s[0], &s[1]);
        }

    auto idx = fwdpp::recycle_mutation_helper(recycling_bin, mutations,
                                              std::move(s), std::move(h), pos,
                                              generation, x);
    lookup.emplace(pos, idx);
    return idx;
}

// For testing, we'll define a single-locus pop type in terms of the above.
using mtype = TwoDmutation;
#define DIPLOID_POPULATION_SIM
#include "common_ind.hpp"

int
main(int argc, char **argv)
{
    TwoDmutation::array_t s{ 0, 1 }, h{ 1, 1 };

    TwoDmutation m(s, h, 0.1, 0, 3);

    if (m.s != s)
        {
            throw std::runtime_error("s not assigned properly!");
        }

    if (m.h != h)
        {
            throw std::runtime_error("h not assigned properly!");
        }

    if (m.pos != 0.1)
        {
            throw std::runtime_error("pos not assigned properly!");
        }

    if (m.neutral)
        {
            throw std::runtime_error("this mutation is not neutral!!");
        }

    if (m.xtra != 3)
        {
            throw std::runtime_error("xtra field set incorrectly");
        }

    std::ostringstream obuffer;
    fwdpp::io::serialize_mutation<TwoDmutation> serializer;
    serializer(obuffer, m);

    std::istringstream ibuffer(obuffer.str());
    fwdpp::io::deserialize_mutation<TwoDmutation> deserializer;
    auto m2 = deserializer(ibuffer);

    if (m != m2)
        {
            throw std::runtime_error("comparison failure");
        }

    // Define a mutation model with no covariance and equal sd in both
    // dimensions

    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, 42);

    diploid_population pop(1000);
    unsigned generation = 0;
    const auto isotropic_mutation =
        [r, &pop, &generation](fwdpp::flagged_mutation_queue &mutation_recbin,
                               diploid_population::mcont_t &mutations) {
            return infsites_TwoDmutation_additive(
                mutation_recbin, mutations, r, pop.mut_lookup, generation, 1.0,
                [r]() { return gsl_rng_uniform(r); }, 0.1, 0.1, 0.0);
        };

    static_assert(
        fwdpp::traits::is_mutation_model<decltype(isotropic_mutation),
                                         diploid_population::mcont_t,
                                         diploid_population::gcont_t>::value,
        "Error: isotropic_mutation is not a valid mutation function");

    auto recbin = fwdpp::make_mut_queue(pop.mcounts);

    auto new_mutation_index = isotropic_mutation(recbin, pop.mutations);

    std::cout << pop.mutations.at(new_mutation_index).pos << ' '
              << pop.mutations.at(new_mutation_index).g << ' '
              << pop.mutations.at(new_mutation_index).neutral << " ["
              << pop.mutations.at(new_mutation_index).s[0] << ','
              << pop.mutations.at(new_mutation_index).s[1] << "] ["
              << pop.mutations.at(new_mutation_index).h[0] << ','
              << pop.mutations.at(new_mutation_index).h[1] << "]\n";

    gsl_rng_free(r);
}
