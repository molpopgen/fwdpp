//  -*- C++ -*-
#ifndef __FWDPP_RECOMBINATION_TCC__
#define __FWDPP_RECOMBINATION_TCC__

#include <vector>
#include <fwdpp/type_traits.hpp>
#include <fwdpp/internal/recombination_common.hpp>
#include <fwdpp/internal/recycling.hpp>
namespace KTfwd
{
    template <typename vec_t, typename gcont_t, typename mcont_t,
              typename queue_t>
    std::size_t
    recombine_gametes(
        const vec_t &pos, gcont_t &gametes, const mcont_t &mutations,
        const std::size_t g1, const std::size_t g2,
        queue_t &gamete_recycling_bin,
        typename gcont_t::value_type::mutation_container &neutral,
        typename gcont_t::value_type::mutation_container &selected)
    {
        assert(g1 < gametes.size());
        assert(g2 < gametes.size());
        assert(std::is_sorted(pos.begin(), pos.end()));
        assert(!pos.empty());
        assert(*(pos.end() - 1) == std::numeric_limits<double>::max());

        // We defer clearing all the way to this point
        // Note that this just adjusts the value of end,
        // and capacity() is not affected, meaning
        // we can insert into these containers relatively
        // efficiently
        neutral.clear();
        selected.clear();
        /*
          Fill neutral and selected by recombining gametes[g1] and gametes[g2]
          at positions contains in pos
        */
        fwdpp_internal::recombine_gametes(pos, g1, g2, gametes, mutations,
                                          neutral, selected);
        return fwdpp_internal::recycle_gamete(gametes, gamete_recycling_bin,
                                              neutral, selected);
    }

    template <typename gcont_t, typename mcont_t, typename recbin_t,
              typename recpol_t>
    std::pair<std::size_t, unsigned>
    recombination(gcont_t &gametes, recbin_t &gamete_recycling_bin,
                  typename gcont_t::value_type::mutation_container &neutral,
                  typename gcont_t::value_type::mutation_container &selected,
                  const recpol_t &rec_pol, const std::size_t g1,
                  const std::size_t g2, const mcont_t &mutations)
    {
        static_assert(
            traits::is_rec_model<recpol_t, typename gcont_t::value_type,
                                 mcont_t>::value,
            "type recpol_t is not a valid recombination policy");
        if (g1 == g2)
            {
                return std::make_pair(g1, 0u);
            }
        auto nm1
            = gametes[g1].mutations.size() + gametes[g1].smutations.size();
        auto nm2
            = gametes[g2].mutations.size() + gametes[g2].smutations.size();
        if ((std::min(nm1, nm2) == 0 && std::max(nm1, nm2) == 1))
            {
                return std::make_pair(g1, 0u);
            }
        /*
          This check is relatively expensive, but we do it for several reasons.
          1. If it returns true, then it was cheaper than a bunch of calls to
          std::upper_bound that would have resulted in nothing being done.
          2. It keeps the number of calls to the RNG the same, meaning that
          the output is the same as previous library versions.
        */
        if (gametes[g1] == gametes[g2])
            return std::make_pair(g1, 0u);
        auto pos = rec_pol(gametes[g1], gametes[g2], mutations);
        if (pos.empty())
            {
                return std::make_pair(g1, 0u);
            }

        assert(pos.back() == std::numeric_limits<double>::max());
        return std::make_pair(recombine_gametes(pos, gametes, mutations, g1,
                                                g2, gamete_recycling_bin,
                                                neutral, selected),
                              unsigned(pos.size() - 1));
    }
}

#endif
