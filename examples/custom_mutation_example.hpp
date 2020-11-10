#pragma once

#include <tuple>
#include <fwdpp/type_traits.hpp>
#include <fwdpp/fundamental_types/mutation_base.hpp>

/*
 * The simplest mutation type, adding just a selection coefficient
 * and dominance to the interface.
*/
struct mutation : public fwdpp::mutation_base
{
    // selection coefficient
    double s;
    // dominance coefficient
    double h;
    mutation(const double &position, const double &sel_coeff,
             const double &dominance = 0.5) noexcept
        : mutation_base(position, (sel_coeff == 0)), s(sel_coeff), h(dominance)
    {
    }

    bool
    operator==(const mutation &rhs) const
    {
        return std::tie(this->s, this->h) == std::tie(rhs.s, rhs.h) && is_equal(rhs);
    }
};

static_assert(fwdpp::traits::is_mutation_v<mutation>,
              "Mutation is not a valid mutation type!");
