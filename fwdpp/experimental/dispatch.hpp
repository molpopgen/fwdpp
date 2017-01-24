#ifndef FWDPP_EXPERIMENTAL_DISPATCH_HPP
#define FWDPP_EXPERIMENTAL_DISPATCH_HPP
/*!
  \file fwdpp/experimental/dispatch.hpp

  \brief SFINAE dispatch of operations involving "custom rules"
*/

#include <type_traits>
#include <functional>

namespace KTfwd
{
    namespace experimental
    {
        template <typename rules_t, typename dipcont_t, typename gcont_t,
                  typename mcont_t, typename fitness_func_t>
        inline auto
        dispatch_w(rules_t &&r, const dipcont_t &diploids, gcont_t &gametes,
                   const mcont_t &mutations, const fitness_func_t &ff)
            -> decltype(r.w(diploids, gametes, mutations, ff))
        {
            r.w(diploids, gametes, mutations, ff);
        }

        template <typename rules_t, typename dipcont_t, typename gcont_t,
                  typename mcont_t, typename fitness_func_t>
        inline auto
        dispatch_w(rules_t &&r, const dipcont_t &diploids, gcont_t &gametes,
                   const mcont_t &mutations, const fitness_func_t &)
            -> decltype(r.w(diploids, gametes, mutations))
        {
            r.w(diploids, gametes, mutations);
        }

        template <typename rules_t, typename diploid_t, typename gcont_t,
                  typename mcont_t, typename fitness_func_t>
        inline auto
        dispatch_update(rules_t &&r, const gsl_rng *rng, diploid_t &offspring,
                        const diploid_t &parent1, const diploid_t &parent2,
                        const gcont_t &gametes, const mcont_t &mutations,
                        const fitness_func_t &)
            -> decltype(r.update(rng, offspring, parent1, parent2, gametes,
                                 mutations))
        {
            r.update(rng, offspring, parent1, parent2, gametes, mutations);
        }

        template <typename rules_t, typename diploid_t, typename gcont_t,
                  typename mcont_t, typename fitness_func_t>
        inline auto
        dispatch_update(rules_t &&r, const gsl_rng *rng, diploid_t &offspring,
                        const diploid_t &parent1, const diploid_t &parent2,
                        const gcont_t &gametes, const mcont_t &mutations,
                        const fitness_func_t &ff)
            -> decltype(r.update(rng, offspring, parent1, parent2, gametes,
                                 mutations, ff))
        {
            r.update(rng, offspring, parent1, parent2, gametes, mutations, ff);
        }
    }
}

#endif
