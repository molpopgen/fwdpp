#ifndef __KTFWD_DEBUG_HPP__
#define __KTFWD_DEBUG_HPP__

#include <algorithm>
#include <numeric>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/type_traits.hpp>

namespace KTfwd
{

    /*! \brief Returns true if the sum of counts in gametes equals twoN, false
      otherwise
      Returns true if the sum of counts in gametes equals twoN, false otherwise
    */
    template <typename gcont_t>
    bool
    check_sum(const gcont_t &gametes, const unsigned twoN)
    {
        static_assert(
            typename traits::is_gamete<typename gcont_t::value_type>::type(),
            "gcont_t::value_type must be a valid gamete type");
        return (std::accumulate(gametes.cbegin(), gametes.cend(), 0u,
                                [](unsigned &__u,
                                   const typename gcont_t::value_type &__g) {
                                    return __u + __g.n;
                                })
                == twoN);
    }

    /*! \brief Returns true if the sum of counts in gametes equals twoN, false
      otherwise
      Returns true if the sum of counts in gametes equals twoN, false otherwise
    */
    template <typename gcont_t>
    bool
    check_sum(const gcont_t *gametes, const unsigned twoN)
    {
        return check_sum(*gametes, twoN);
    }

    template <typename gamete_t, typename mcont_t>
    bool
    gamete_is_sorted_n(const gamete_t &g, const mcont_t &m)
    /*!
      \brief Check that neutral mutation keys are sorted according to mutation
      position
    */
    {
        return std::is_sorted(g.mutations.begin(), g.mutations.end(),
                              [&m](const size_t i, const size_t j) {
                                  return m[i].pos <= m[j].pos;
                              });
    }

    template <typename gamete_t, typename mcont_t>
    bool
    gamete_is_sorted_s(const gamete_t &g, const mcont_t &m)
    /*!
      \brief Check that selected mutation keys are sorted according to mutation
      position
    */
    {
        return std::is_sorted(g.smutations.begin(), g.smutations.end(),
                              [&m](const size_t i, const size_t j) {
                                  return m[i].pos <= m[j].pos;
                              });
    }

    template <typename gamete_t, typename mcont_t>
    bool
    gamete_data_sane(const gamete_t &g, const mcont_t &mutations,
                     const std::vector<uint_t> &mutcounts)
    /*
      \brief Check that "neutral" and "non-neutral" mutations are where we
      expect them to be.
     */
    {
        for (const auto &i : g.mutations)
            {
                if (!mutcounts[i])
                    return false;
                if (!mutations[i].neutral)
                    return false;
                if (!(g.n <= mutcounts[i]))
                    return false;
            }
        for (const auto &i : g.smutations)
            {
                if (!mutcounts[i])
                    return false;
                if (mutations[i].neutral)
                    return false;
                if (!(g.n <= mutcounts[i]))
                    return false;
            }
        return true;
    }

    template <typename dipcont_t, typename gcont_t, typename mcont_t>
    bool
    popdata_sane(const dipcont_t &diploids, const gcont_t &gametes,
                 const mcont_t &mutations,
                 const std::vector<uint_t> &mutcounts)
    /*
      \brief Check that all diploids refer to extant gametes with sane data.
     */
    {
        for (const auto &d : diploids)
            {
                if (!gametes[d.first].n)
                    return false;
                if (!gametes[d.second].n)
                    return false;
                if (!gamete_data_sane(gametes[d.first], mutations, mutcounts))
                    return false;
                if (!gamete_data_sane(gametes[d.second], mutations, mutcounts))
                    return false;
            }
        return true;
    }

    template <typename dipcont_t, typename gcont_t, typename mcont_t>
    bool
    popdata_sane_multilocus(const dipcont_t &diploids, const gcont_t &gametes,
                            const mcont_t &mutations,
                            const std::vector<uint_t> &mutcounts)
    /*
      \brief Check that all diploids refer to extant gametes with sane data.
      \note This is for the case where diploids are equivalent to
      vector<pair<key,key>>
     */
    {
        for (const auto &d : diploids)
            {
                for (const auto &locus : d)
                    {
                        if (!gametes[locus.first].n)
                            return false;
                        if (!gametes[locus.second].n)
                            return false;
                        if (!gamete_data_sane(gametes[locus.first], mutations,
                                              mutcounts))
                            return false;
                        if (!gamete_data_sane(gametes[locus.second], mutations,
                                              mutcounts))
                            return false;
                    }
            }
        return true;
    }
}

#endif
