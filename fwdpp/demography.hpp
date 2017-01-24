/*!
  \file demography.hpp

  \brief Low-level functions for demographic events
*/

/*!
  \defgroup demography

  \brief Functions relating to demographic events.
*/

#ifndef FWDPP_DEMOGRAPHY_HPP
#define FWDPP_DEMOGRAPHY_HPP

#include <cassert>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <set>
#include <gsl/gsl_rng.h>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/internal/sample_diploid_helpers.hpp>
#include <fwdpp/internal/demography_details.hpp>

namespace KTfwd
{
    /*! \brief Copy a deme

      Make an exact copy of the i-th deme.  The new deme is placed at the end
      of \c diploids.

      \param mutations Container of mutations
      \param mcounts Container of mutation counts
      \param gametes Container of gametes
      \param diploids Container of diploids in a metapopulation (vector of
      vectors of diploids).
      \param i The index of the deme to copy.

      \return An error value -1 if i >= diploids.size(), 0 otherwise (meaning
      copy successful).

      \note Updates gamete and mutation counts

      \ingroup demography
    */
    template <typename mcont_t, typename mcount_t, typename gcont_t,
              typename vdipvector_t>
    int
    copy_deme(const mcont_t &mutations, mcount_t &mcounts, gcont_t &gametes,
              vdipvector_t &diploids, const size_t i)
    {
        if (i >= diploids.size())
            return -1;
        diploids.emplace_back(diploids[i]);
        // We've added a deme, so we need to go through it, and
        // update gamete counts accordingly
        for (auto &dip : diploids[diploids.size() - 1])
            {
                gametes[dip.first].n++;
                gametes[dip.second].n++;
            }
        fwdpp_internal::process_gametes(gametes, mutations, mcounts);
        return 0;
    }

    /*! \brief Merge demes

      Merge two demes

      \param diploids Container of diploids in a metapopulation (vector of
      vectors of diploids).
      \param i One deme to merge
      \param j The other deme to merge

      Deme max(i,j) is appended to the end of deme min(i,j).  Deme max(i,j) is
      then removed from diploids.

      \return Error value -1 if i or j >= diploids.size(), Error value 1 if
      i==j, 0 otherwise (meaning merge successful).

      \ingroup demography
    */
    template <typename vdipvector_t>
    int
    merge_demes(vdipvector_t &diploids, size_t i, size_t j)
    {
        if (i >= diploids.size() || j >= diploids.size())
            return -1;
        if (i == j)
            return 1;
        if (i > j)
            std::swap(i, j);
        std::move(diploids[j].begin(), diploids[j].end(),
                  std::back_inserter(diploids[i]));
        diploids.erase(diploids.begin() + j);
        return 0;
    }

    /*! \brief Delete a deme

      Delete a deme from a metapopulation

      \param mutations Container of mutations
      \param mcounts Container of mutation counts
      \param gametes Container of gametes
      \param diploids Container of diploids in a metapopulation (vector of
      vectors of diploids).
      \param i The index of the deme to delete.

      \return Error value -1 if i >= diploids.size(), 0 otherwise (meaning
      successful delete).

      \note Gamete and mutation counts are updated.

      \ingroup demography
    */
    template <typename mcont_t, typename mcount_t, typename gcont_t,
              typename vdipvector_t>
    int
    remove_deme(const mcont_t &mutations, mcount_t &mcounts, gcont_t &gametes,
                vdipvector_t &diploids, const size_t i)
    {
        if (i >= diploids.size())
            return -1;
        for (auto &dip : diploids[i]) // update gamete counts
            {
                gametes[dip.first].n--;
                gametes[dip.second].n--;
            }
        diploids.erase(diploids.begin() + i);
        KTfwd::fwdpp_internal::process_gametes(
            gametes, mutations, mcounts); // update mutation counts
        return 0;
    }

    /*! \brief Swap demes

      Swap two demes

      \param diploids Container of diploids in a metapopulation (vector of
      vectors of diploids).
      \param i One deme to swap
      \param j The other deme to swap

      \return Error value -1 if i or j >= diploids.size(), 1 if i == j, 0
      otherwise (meaning successful swap).

      \note Implemented as std::swap(diploids[i],diploids[j]).

      \ingroup demography
     */
    template <typename vdipvector_t>
    int
    swap_demes(vdipvector_t &diploids, const size_t i, const size_t j)
    {
        if (i >= diploids.size() || j >= diploids.size())
            return -1;
        if (i == j)
            return 1;
        std::swap(diploids[i], diploids[j]);
        return 0;
    }

    /*! \brief Split deme into two.

      Take deme \c i and create a new deme of size \c N_new from it.  Deme \c
      i's
      size is reduced \c N_new individuals.

      \param r A pointer to a gsl_rng object

      \param mutations Container of mutations
      \param mcounts Container of mutation counts
      \param gametes Container of gametes
      \param diploids Container of diploids in a metapopulation (vector of
      vectors of diploids).
      \param i The index of the deme to split.
      \param N_new The population size of the new deme.
      \param replacement Sample individuals with replacement or not.

      The behavior of this function depends on the value of \c replacement.  If
      \c replacement is \c false, then \c N_new
      individuals are selected at random (w/o replacement), removed from deme
      \c i, and added into the new deme.  For this case,
      N_new must also be <= diploids[i].size()

      When \c replacement is \c true, then deme \c i is first copied to a
      location outside of \c diploids.  Then, \c N_new individuals are
      randomly chosen (with replacement) from the copy of deme \c i and added
      into the new deme.  Likewise, \c diploids[i].size()-N_new
      individuals are randomly chosen (with replacement) from the copy of deme
      \c i and used to re-populate deme \c i.

      \return Error code -1 if i>=diploids.size(), 1 if (!replacement&&(N_new
      >= diploids[i].size())), 0 otherwise (meaning successful split).

      \note Gamete and mutation counts are updated only when sampling with
      replacement.  Counts cannot change if sampling without replacement.

      \ingroup demography
    */
    template <typename mcont_t, typename mcount_t, typename gcont_t,
              typename vdipvector_t>
    int
    split_deme(const gsl_rng *r, const mcont_t &mutations, mcount_t &mcounts,
               gcont_t &gametes, vdipvector_t &diploids, const size_t i,
               const uint_t N_new, const bool replacement = false)
    {
        if (i >= diploids.size())
            return -1;
        if (!replacement && (N_new >= diploids[i].size()))
            return 1;
        return fwdpp_internal::split_deme_details(
            r, mutations, mcounts, gametes, diploids, i, N_new, replacement);
    }

    /*!
      \brief Create admixed population

       Create an admixed population with a specified amount of ancestry from
      each parental deme.

       \param r A pointer to a gsl_rng object
       \param mutations Container of mutations
       \param mcounts Container of mutation counts
       \param gametes Container of gametes
       \param diploids Container of diploids in a metapopulation (vector of
      vectors of diploids).
       \param i The index of one parental deme
       \param j The index of the other parental deme
       \param pi The fraction of ancestors in the admixed deme coming from deme
      \c i
       \param N_new The size of the new, admixed deme
       \param replacement Whether or not to sample parents with replacement
      from their respective demes.

       The behavior of this function depends on the value of \c replacement.
      When \c replacement is \c false, the size of the
       new deme cannot be so large that it requires all of demes \c i and/or \c
      j.  When \c replacement is \c true, there is no limit
       on the value of \c N_new.

       \return Error value -1 if \c i or \c j >= \c diploids.size(), -1 if \c
      replacement is true and the new deme size is too large,
       1 if \c pi < 0. or \c pi >= 1., 0 otherwise (meaning successful
      admixture).

       \ingroup demography
    */
    template <typename mcont_t, typename mcount_t, typename gcont_t,
              typename vdipvector_t>
    int
    admix_demes(const gsl_rng *r, const mcont_t &mutations, mcount_t &mcounts,
                gcont_t &gametes, vdipvector_t &diploids, const size_t i,
                const size_t j, const double pi, const uint_t N_new,
                const bool replacement = false)
    {
        if (i >= diploids.size() || j >= diploids.size())
            return -1;
        if (pi < 0. || pi >= 1.)
            return 1;

        uint_t N_from_i = std::round(pi * double(N_new)),
               N_from_j = N_new - N_from_i;

        if (!replacement) // check for logic errors in input
            {
                if (N_from_i >= diploids[i].size())
                    return 1;
                if (N_from_j >= diploids[j].size())
                    return 1;
            }

        // add a new deme
        diploids.emplace_back(typename vdipvector_t::value_type());
        // get reference to new deme
        auto &new_deme = diploids[diploids.size() - 1];
        new_deme.reserve(N_new);
        const auto &parental_deme_i = diploids[i];

        // get list of individuals from deme i
        auto indlist = fwdpp_internal::sample_individuals(
            r, diploids[i].size(), N_from_i, replacement);

        // add individuals from deme i into new deme
        for (const auto &ind : indlist)
            new_deme.push_back(parental_deme_i[ind]);

        // get list of individuals from deme j
        indlist = fwdpp_internal::sample_individuals(r, diploids[j].size(),
                                                     N_from_j, replacement);

        // update reference to parental deme
        const auto &parental_deme_j = diploids[j];

        // add individuals from deme i into new deme
        for (const auto &ind : indlist)
            new_deme.push_back(parental_deme_j[ind]);

        // update gamete counts due to new deme
        for (const auto &dip : new_deme)
            {
                gametes[dip.first].n++;
                gametes[dip.second].n++;
            }

        fwdpp_internal::process_gametes(gametes, mutations, mcounts);
        return 0;
    }
}

#endif
