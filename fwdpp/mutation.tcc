//  -*- C++ -*-
#ifndef _FWDPP_MUTATION_TCC_
#define _FWDPP_MUTATION_TCC_

#include <type_traits>
#include <numeric>
#include <gsl/gsl_randist.h>
#include <fwdpp/internal/mutation_internal.hpp>

namespace KTfwd
{
    template <typename queue_type, typename queue_type2,
              typename mutation_model, typename gamete_insertion_policy,
              typename gcont_t, typename mcont_t>
    std::size_t
    mutate_gamete_recycle(queue_type &recycling_bin,
                          queue_type2 &gamete_recycling_bin, const gsl_rng *r,
                          const double &mu, gcont_t &gametes,
                          mcont_t &mutations, const std::size_t g,
                          const mutation_model &mmodel,
                          const gamete_insertion_policy &gpolicy)
    {
        static_assert(traits::is_mutation_model<mutation_model, mcont_t,
                                                   gcont_t>::value,
                      "error: type mutation_model is not a dispatchable "
                      "mutation model type!");
        assert(g < gametes.size());
        assert(gamete_is_sorted_n(gametes[g], mutations));
        assert(gamete_is_sorted_s(gametes[g], mutations));
        unsigned nm = gsl_ran_poisson(r, mu);
        if (!nm)
            return g;

        assert(gametes[g].n);
        gametes[g].n--;
        // Recycle an already-allocated gamete, if possible
        if (!gamete_recycling_bin.empty())
            {
                auto idx = gamete_recycling_bin.front();
                gamete_recycling_bin.pop();
                assert(idx != g);
                assert(!gametes[idx].n);
                gametes[idx].mutations = gametes[g].mutations;
                gametes[idx].smutations = gametes[g].smutations;
                gametes[idx].n = 1;
                fwdpp_internal::add_N_mutations_recycle(
                    recycling_bin, mmodel, nm, mutations, gametes[idx]);
                return idx;
            }
        // If not, create a new gamete that we'll add to the gamete container
        typename gcont_t::value_type ng(1, gametes[g].mutations,
                                        gametes[g].smutations);
        fwdpp_internal::add_N_mutations_recycle(recycling_bin, mmodel, nm,
                                                mutations, ng);
        return gpolicy(std::move(ng), gametes);
    }
}
#endif /* _FWDPP_MUTATION_TCC_ */
