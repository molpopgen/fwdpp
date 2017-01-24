/*!
  \file fwdpp/sugar/add_mutation.hpp

  \brief Definition of add_mutation and add_mutations
*/

#ifndef FWDPP_SUGAR_ADD_MUTATION_HPP
#define FWDPP_SUGAR_ADD_MUTATION_HPP

#include <exception>
#include <type_traits>
#include <algorithm>
#include <cassert>
#include <vector>
#include <unordered_map>
#include <fwdpp/internal/recycling.hpp>
#include <fwdpp/sugar/poptypes/tags.hpp>

namespace KTfwd
{
    namespace sugar
    {
        template <typename mcont_t, typename mcounts_t>
        std::size_t
        get_mut_index(mcont_t &mutations, mcounts_t &mcounts,
                      typename mcont_t::value_type &new_mutation)
        /*!
          \brief Helper function for implementing KTfwd::add_mutation

          If \a new_mutation already exists in \a mutations, its location is
          returned.

          Otherwise, this function puts the new mutation into the mutations
          container and
          updates other objects as needed.

          \return A key telling where \a new_mutation exists within \a
          mutations.

          \note This function requires that operator== be defined for the
          mutation type.

          \note In the event that \a new_mutation is added to \a mutations, the
          count for that
          mutations will be zero.

          \ingroup sugar
        */
        {
            // If the mutation already exists, let's not create it...
            auto mitr = std::find(std::begin(mutations), std::end(mutations),
                                  new_mutation);
            if (mitr != mutations.end())
                return std::size_t(std::distance(std::begin(mutations), mitr));

            /*
              If we are here, then new_mutation does not currently exist, so we
              have to add it.

              We do this using fwdpp's recycling method.  But, we don't need to
              create recycling bin,
              as we're just looking to add one mutation.  Instead, we can
              simply search for an extinct
              mutation.
            */
            // find an extinct mutation, if one exists
            auto extinct_mut
                = std::find(std::begin(mcounts), std::end(mcounts), 0);
            std::size_t mindex = mcounts.size(); // this will be the correct
                                                 // value if we use the else
                                                 // block below
            if (extinct_mut != std::end(mcounts)) // then we can recycle
                {
                    auto dist
                        = std::distance(std::begin(mcounts), extinct_mut);
                    mutations[dist] = std::move(
                        new_mutation);          // move new mutation into place
                    mindex = std::size_t(dist); // update our value
                }
            else
                {
                    // cannot recycle, so add it to end
                    mutations.emplace_back(std::move(new_mutation));
                    mcounts.push_back(0); // Add a place for this variant
                }
            return mindex;
        }

        template <typename poptype, typename map_t>
        void
        add_mutation_details(poptype &p,
                             const std::vector<std::size_t> &mindexes,
                             const map_t &gams)
        {
            auto gam_recycling_bin
                = fwdpp_internal::make_gamete_queue(p.gametes);
            // Function object for calls to upper bound
            auto inserter
                = [&p](const double &__value, const std::size_t __mut) noexcept
            {
                assert(__mut < p.mutations.size());
                return __value < p.mutations[__mut].pos;
            };

            for (const auto &gi : gams)
                {
                    auto n = p.gametes[gi.first].mutations;
                    auto s = p.gametes[gi.first].smutations;
                    for (auto mindex : mindexes)
                        {
                            auto pos = p.mutations[mindex].pos;

                            if (p.mutations[mindex].neutral)
                                {
                                    n.insert(std::upper_bound(n.begin(),
                                                              n.end(), pos,
                                                              inserter),
                                             mindex);
                                }
                            else
                                {
                                    s.insert(std::upper_bound(s.begin(),
                                                              s.end(), pos,
                                                              inserter),
                                             mindex);
                                }
                            // update mutation count
                            p.mcounts[mindex] +=
                                typename poptype::mcount_t::value_type(
                                    gi.second.size());
                        }
                    // get new gamete
                    auto new_gamete_key = fwdpp_internal::recycle_gamete(
                        p.gametes, gam_recycling_bin, n, s);
                    // update gamete count
                    p.gametes[gi.first].n
                        -= decltype(p.gametes[gi.first].n)(gi.second.size());
                    p.gametes[new_gamete_key].n += decltype(
                        p.gametes[new_gamete_key].n)(gi.second.size());

                    // This updates every diploid to have a key to
                    // p.gametes[new_gamete_key]
                    for (auto i : gi.second)
                        {
                            *i = new_gamete_key;
                        }
                }
        }

        template <typename poptype, typename = std::enable_if<std::is_same<
                                        typename poptype::popmodel_t,
                                        KTfwd::sugar::SINGLEPOP_TAG>::value>>
        std::unordered_map<std::size_t,
                           std::vector<
                               typename poptype::diploid_t::first_type *>>
        collect_gametes(poptype &p, const std::vector<std::size_t> &indlist,
                        const std::vector<short> &clist)
        /*!
          \brief Helper function for add_mutation and add_mutations.

          Collects together all the gametes that will be updated, to minimize
          the number of
          new gametes created.
        */
        {
            std::unordered_map<std::size_t,
                               std::vector<
                                   typename poptype::diploid_t::first_type *>>
                gams;
            for (std::size_t i = 0; i < indlist.size(); ++i)
                {
                    if (clist[i] == 0 || clist[i] == 2)
                        {
                            gams[p.diploids[indlist[i]].first].push_back(
                                &p.diploids[indlist[i]].first);
                        }
                    if (clist[i] > 0)
                        {
                            gams[p.diploids[indlist[i]].second].push_back(
                                &p.diploids[indlist[i]].second);
                        }
                }
            return gams;
        }

        template <typename poptype, typename = std::enable_if<std::is_same<
                                        typename poptype::popmodel_t,
                                        KTfwd::sugar::MULTILOCPOP_TAG>::value>>
        std::unordered_map<std::size_t,
                           std::vector<typename poptype::diploid_t::
                                           value_type::first_type *>>
        collect_gametes(poptype &p, const std::size_t locus,
                        const std::vector<std::size_t> &indlist,
                        const std::vector<short> &clist)
        /*!
          \brief Helper function for add_mutation and add_mutations.

          Collects together all the gametes that will be updated, to minimize
          the number of
          new gametes created.
        */
        {
            std::unordered_map<std::size_t,
                               std::vector<typename poptype::diploid_t::
                                               value_type::first_type *>>
                gams;
            for (std::size_t ind = 0; ind < indlist.size(); ++ind)
                {
                    if (clist[ind] == 0 || clist[ind] == 2)
                        {
                            gams[p.diploids[indlist[ind]][locus].first]
                                .push_back(
                                    &p.diploids[indlist[ind]][locus].first);
                        }
                    if (clist[ind] > 0)
                        {
                            gams[p.diploids[indlist[ind]][locus].second]
                                .push_back(
                                    &p.diploids[indlist[ind]][locus].second);
                        }
                }
            return gams;
        }

        template <typename poptype, typename = std::enable_if<std::is_same<
                                        typename poptype::popmodel_t,
                                        KTfwd::sugar::METAPOP_TAG>::value>>
        std::unordered_map<std::size_t,
                           std::vector<
                               typename poptype::diploid_t::first_type *>>
        collect_gametes(poptype &p, const std::vector<std::size_t> demes,
                        const std::vector<std::vector<std::size_t>> &indlist,
                        const std::vector<std::vector<short>> &clist)
        /*!
          \brief Helper function for add_mutation and add_mutations.

          Collects together all the gametes that will be updated, to minimize
          the number of
          new gametes created.
        */
        {
            // collect gametes that are to be updated
            std::unordered_map<std::size_t,
                               std::vector<
                                   typename poptype::diploid_t::first_type *>>
                gams;
            for (std::size_t deme = 0; deme < demes.size(); ++deme)
                {
                    for (std::size_t ind = 0; ind < indlist[deme].size();
                         ++ind)
                        {
                            auto c = clist[deme][ind];
                            if (c == 0 || c == 2)
                                {
                                    gams[p.diploids[demes[deme]]
                                                   [indlist[deme][ind]]
                                                       .first]
                                        .push_back(
                                            &p.diploids[demes[deme]]
                                                       [indlist[deme][ind]]
                                                           .first);
                                }
                            if (c > 0)
                                {
                                    gams[p.diploids[demes[deme]]
                                                   [indlist[deme][ind]]
                                                       .second]
                                        .push_back(
                                            &p.diploids[demes[deme]]
                                                       [indlist[deme][ind]]
                                                           .second);
                                }
                        }
                }
            return gams;
        }
    }

    template <typename poptype, class... Args>
    std::size_t
    add_mutation(poptype &p, const std::vector<std::size_t> &indlist,
                 const std::vector<short> &clist, Args &&... args)
    /*!
      \brief Add a mutation into a population at a given frequency.

      \param p A single deme object.
      \param indlist A list of indexes of diploids into which to add the new
      mutations.
      \param clist A list of gametes.  See below.
      \param args Values required to cosnstruct a new mutation.  See below.

      \return The key referring to the location of the new mutation in the
      population

      Some notes:

      clist.size() must equal indlist.size()

      Values in \a clist must be 0, 1, or 2. These values mean to add the
      mutation to the first,
      second, or both gametes, resepectively, of each diploid in \a indlist.

      Note that \a args can take on a few different forms.  First, it can be a
      raw set of values
      used to construct a new mutation.  Or, it can be an object of correct
      mutation type.  Or, it can be
      any type from which the correct mutation type can be constructed.  The
      last two cases require
      that the mutation type have the appropriate constructors defined.

      See the unit test file unit/test_sugar_add_mutation.cc for example of
      use.

      \ingroup sugar
    */
    {
        static_assert(std::is_same<typename poptype::popmodel_t,
                                   KTfwd::sugar::SINGLEPOP_TAG>::value,
                      "poptype must be a single-deme object type");

        // Before we go deep into creating objects, let's do some checks
        for (const auto &i : indlist)
            {
                if (i >= p.diploids.size())
                    throw std::out_of_range(
                        "indlist contains elements > p.diploids.size()");
            }
        for (const auto &c : clist)
            {
                if (c < 0 || c > 2)
                    throw std::runtime_error(
                        "clist contains elements < 0 and/or > 2");
            }
        if (indlist.size() != clist.size())
            throw std::runtime_error("indlist and clist must be same length");

        // create a new mutation
        typename poptype::mcont_t::value_type new_mutant(
            std::forward<Args>(args)...);
        auto mindex = sugar::get_mut_index(p.mutations, p.mcounts, new_mutant);
        auto gams = collect_gametes(p, indlist, clist);
        sugar::add_mutation_details(p, { mindex }, gams);
        return mindex;
    }

    template <typename metapoptype, class... Args>
    std::size_t
    add_mutation(metapoptype &p, const std::vector<std::size_t> demes,
                 const std::vector<std::vector<std::size_t>> &indlist,
                 const std::vector<std::vector<short>> &clist, Args &&... args)
    /*!
      \brief Add a mutation into a deme from a population at a given frequency
      in a specific set of
      demes.

      \param p A population object. Meta- or multi-locus.
      \param demes Vector of indices of the demes in which to add mutation.
      \param indlist A list of indexes of diploids into which to add the new
      mutations in each deme.
      \param clist A list of gametes. See below.
      \param args Values required to cosnstruct a new mutation.  See below.

      \return The key referring to the location of the new mutation in the
      population

      Some notes:

      clist.size() must equal indlist.size() must equal demes.size()

      Further, each element of clist must equal each element of indlist in
      size.

      Values in \a clist must be 0, 1, or 2. These values mean to add the
      mutation to the first,
      second, or both gametes, resepectively, of each diploid in \a indlist.

      Note that \a args can take on a few different forms.  First, it can be a
      raw set of values
      used to construct a new mutation.  Or, it can be an object of correct
      mutation type.  Or, it can be
      any type from which the correct mutation type can be constructed.  The
      last two cases require
      that the mutation type have the appropriate constructors defined.

      See the unit test file unit/test_sugar_add_mutation.cc for example of
      use.

      \ingroup sugar
    */
    {
        for (const auto &c : clist)
            {
                for (auto ci : c)
                    if (ci < 0 || ci > 2)
                        throw std::runtime_error(
                            "clist contains elements < 0 and/or > 2");
            }
        if (indlist.size() != clist.size())
            throw std::runtime_error("indlist and clist must be same length");
        for (std::size_t i = 0; i < indlist.size(); ++i)
            {
                if (indlist[i].size() != clist[i].size())
                    throw std::runtime_error(
                        "indlist[i] must equal clist[i] for all i");
            }

        for (std::size_t i = 0; i < demes.size(); ++i)
            {
                if (i > p.diploids.size())
                    throw std::out_of_range("deme index out of range");
                auto d = demes[i];
                for (const auto di : indlist[i])
                    {
                        if (di > p.diploids[d].size())
                            throw std::out_of_range(
                                "individual index out of range");
                    }
            }
        // create a new mutation
        typename metapoptype::mcont_t::value_type new_mutant(
            std::forward<Args>(args)...);
        auto mindex = sugar::get_mut_index(p.mutations, p.mcounts, new_mutant);
        auto gams = sugar::collect_gametes(p, demes, indlist, clist);
        sugar::add_mutation_details(p, { mindex }, gams);
        return mindex;
    }

    template <typename multiloc_poptype, class... Args>
    std::size_t
    add_mutation(multiloc_poptype &p, const std::size_t locus,
                 const std::vector<std::size_t> &indlist,
                 const std::vector<short> &clist, Args &&... args)
    /*!
      \brief Add a mutation into a population at a given frequency at in a
      specific locus.

      \param p A population object. Meta- or multi-locus.
      \param locus Index of the locus in which to add mutation.
      \param indlist A list of indexes of diploids into which to add the new
      mutations.
      \param clist A list of gametes. See below.
      \param args Values required to cosnstruct a new mutation.  See below.


      \return Nothing (void)    \return The key referring to the location of
      the new mutation in the population

      Some notes:

      clist.size() must equal indlist.size()

      Values in \a clist must be 0, 1, or 2. These values mean to add the
      mutation to the first,
      second, or both gametes, resepectively, of each diploid in \a indlist.

      Note that \a args can take on a few different forms.  First, it can be a
      raw set of values
      used to construct a new mutation.  Or, it can be an object of correct
      mutation type.  Or, it can be
      any type from which the correct mutation type can be constructed.  The
      last two cases require
      that the mutation type have the appropriate constructors defined.

      See the unit test file unit/test_sugar_add_mutation.cc for example of
      use.

      \ingroup sugar
    */
    {
        for (const auto &c : clist)
            {
                if (c < 0 || c > 2)
                    throw std::runtime_error(
                        "clist contains elements < 0 and/or > 2");
            }
        if (indlist.size() != clist.size())
            throw std::runtime_error("indlist and clist must be same length");

        if (locus >= p.diploids[0].size())
            throw std::out_of_range("locus index out of range");
        for (const auto &i : indlist)
            {
                if (i >= p.diploids.size())
                    throw std::out_of_range(
                        "indlist contains elements > p.diploids.size()");
            }

        // create a new mutation
        typename multiloc_poptype::mcont_t::value_type new_mutant(
            std::forward<Args>(args)...);
        auto mindex = sugar::get_mut_index(p.mutations, p.mcounts, new_mutant);
        auto gams = sugar::collect_gametes(p, locus, indlist, clist);
        sugar::add_mutation_details(p, { mindex }, gams);
        return mindex;
    }

    template <typename poptype>
    void
    add_mutations(poptype &p, const std::vector<std::size_t> &indlist,
                  const std::vector<short> &clist,
                  const std::vector<std::size_t> &mutation_indexes)
    /*!
      \brief Add a set of mutations into a set of individuals in a population.

      \param p A population object. Meta- or multi-locus.
      \param indlist A list of indexes of diploids into which to add the new
      mutations.
      \param clist A list of gametes. See below.
      \param mutation_indexes The set of mutations to add.  See below.

      \return Nothing (void)

      Some notes:

      clist.size() must equal indlist.size()

      Values in \a clist must be 0, 1, or 2. These values mean to add the
      mutation to the first,
      second, or both gametes, resepectively, of each diploid in \a indlist.

      \note \a mutation_indexes refers to the locations of mutations found in
      \a p.mutations.

      \note For each element, i, in \a mutation_indexes, \a p.mcounts[i] should
      be zero before
      entering this function.

      \note \a p.mut_lookup is NOT updated by this function.

      See the unit test file unit/test_sugar_add_mutation.cc for example of
      use.

      \ingroup sugar
    */
    {
        static_assert(std::is_same<typename poptype::popmodel_t,
                                   KTfwd::sugar::SINGLEPOP_TAG>::value,
                      "poptype must be a single-deme object type");

        // Before we go deep into creating objects, let's do some checks
        for (const auto &i : indlist)
            {
                if (i >= p.diploids.size())
                    throw std::out_of_range(
                        "indlist contains elements > p.diploids.size()");
            }
        for (const auto &c : clist)
            {
                if (c < 0 || c > 2)
                    throw std::runtime_error(
                        "clist contains elements < 0 and/or > 2");
            }
        if (indlist.size() != clist.size())
            throw std::runtime_error("indlist and clist must be same length");
        for (auto mi : mutation_indexes)
            {
                if (mi >= p.mutations.size())
                    throw std::runtime_error("mutation key out of range");
            }
        if (p.mcounts.size() != p.mutations.size())
            throw std::runtime_error("p.mcounts.size() != p.mutations.size()");
        auto gams = sugar::collect_gametes(p, indlist, clist);
        sugar::add_mutation_details(p, mutation_indexes, gams);
    }

    template <typename metapoptype>
    void
    add_mutations(metapoptype &p, const std::vector<std::size_t> demes,
                  const std::vector<std::vector<std::size_t>> &indlist,
                  const std::vector<std::vector<short>> &clist,
                  const std::vector<std::size_t> &mutation_indexes)
    /*!
      \brief Add a set of mutation into a deme from a population at a given
      frequency in a specific set of
      demes.

      \param p A population object. Meta- or multi-locus.
      \param demes Vector of indices of the demes in which to add mutation.
      \param indlist A list of indexes of diploids into which to add the new
      mutations in each deme.
      \param clist A list of gametes. See below.
      \param mutation_indexes Keys to mutations in p.mutations

      \return Nothing (void)

      Some notes:

      clist.size() must equal indlist.size()

      Values in \a clist must be 0, 1, or 2. These values mean to add the
      mutation to the first,
      second, or both gametes, resepectively, of each diploid in \a indlist.

      \note \a mutation_indexes refers to the locations of mutations found in
      \a p.mutations.

      \note For each element, i, in \a mutation_indexes, \a p.mcounts[i] should
      be zero before
      entering this function.

      \note \a p.mut_lookup is NOT updated by this function.

      See the unit test file unit/test_sugar_add_mutation.cc for example of
      use.

      \ingroup sugar
    */
    {
        for (const auto &c : clist)
            {
                for (auto ci : c)
                    if (ci < 0 || ci > 2)
                        throw std::runtime_error(
                            "clist contains elements < 0 and/or > 2");
            }
        if (indlist.size() != clist.size())
            throw std::runtime_error("indlist and clist must be same length");
        for (std::size_t i = 0; i < indlist.size(); ++i)
            {
                if (indlist[i].size() != clist[i].size())
                    throw std::runtime_error(
                        "indlist[i] must equal clist[i] for all i");
            }

        for (std::size_t i = 0; i < demes.size(); ++i)
            {
                if (i > p.diploids.size())
                    throw std::out_of_range("deme index out of range");
                auto d = demes[i];
                for (const auto di : indlist[i])
                    {
                        if (di > p.diploids[d].size())
                            throw std::out_of_range(
                                "individual index out of range");
                    }
            }
        for (auto mi : mutation_indexes)
            {
                if (mi >= p.mutations.size())
                    throw std::runtime_error("mutation key out of range");
            }
        if (p.mcounts.size() != p.mutations.size())
            throw std::runtime_error("p.mcounts.size() != p.mutations.size()");
        auto gams = sugar::collect_gametes(p, demes, indlist, clist);
        sugar::add_mutation_details(p, mutation_indexes, gams);
    }

    template <typename multiloc_poptype>
    void
    add_mutations(multiloc_poptype &p, const std::size_t locus,
                  const std::vector<std::size_t> &indlist,
                  const std::vector<short> &clist,
                  const std::vector<std::size_t> &mutation_indexes)
    /*!
      \brief Add a set of mutations into a given locus of a multi-locus
      simulation.

      \param p A population object. Meta- or multi-locus.
      \param locus Index of the locus in which to add mutation.
      \param indlist A list of indexes of diploids into which to add the new
      mutations.
      \param clist A list of gametes. See below.
      \param mutaton_indexes Keys to mutations in p.mutations.

      \return Nothing (void)

      Some notes:

      clist.size() must equal indlist.size()

      Values in \a clist must be 0, 1, or 2. These values mean to add the
      mutation to the first,
      second, or both gametes, resepectively, of each diploid in \a indlist.

      \note \a mutation_indexes refers to the locations of mutations found in
      \a p.mutations.

      \note For each element, i, in \a mutation_indexes, \a p.mcounts[i] should
      be zero before
      entering this function.

      \note \a p.mut_lookup is NOT updated by this function.

      See the unit test file unit/test_sugar_add_mutation.cc for example of
      use.

      \ingroup sugar
    */
    {
        for (const auto &c : clist)
            {
                if (c < 0 || c > 2)
                    throw std::runtime_error(
                        "clist contains elements < 0 and/or > 2");
            }
        if (indlist.size() != clist.size())
            throw std::runtime_error("indlist and clist must be same length");

        if (locus >= p.diploids[0].size())
            throw std::out_of_range("locus index out of range");
        for (const auto &i : indlist)
            {
                if (i >= p.diploids.size())
                    throw std::out_of_range(
                        "indlist contains elements > p.diploids.size()");
            }
        for (auto mi : mutation_indexes)
            {
                if (mi >= p.mutations.size())
                    throw std::runtime_error("mutation key out of range");
            }
        if (p.mcounts.size() != p.mutations.size())
            throw std::runtime_error("p.mcounts.size() != p.mutations.size()");
        auto gams = sugar::collect_gametes(p, locus, indlist, clist);
        sugar::add_mutation_details(p, mutation_indexes, gams);
    }
}

#endif
