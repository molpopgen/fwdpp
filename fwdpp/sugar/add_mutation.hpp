/*!
  \file fwdpp/sugar/add_mutation.hpp

  \brief Definition of add_mutation and add_mutations
*/

#ifndef FWDPP_SUGAR_ADD_MUTATION_HPP
#define FWDPP_SUGAR_ADD_MUTATION_HPP

#include <exception>
#include <type_traits>
#include <algorithm>
#include <stdexcept>
#include <vector>
#include <unordered_map>
#include <fwdpp/simfunctions/recycling.hpp>
#include <fwdpp/poptypes/tags.hpp>
#include <fwdpp/debug.hpp>

namespace fwdpp
{

    namespace sugar
    {
        template <typename MutationContainerType, typename mcounts_t>
        std::size_t
        get_mut_index(MutationContainerType &mutations, mcounts_t &mcounts,
                      typename MutationContainerType::value_type &new_mutation)
        /*!
          \brief Helper function for implementing fwdpp::add_mutation

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
            auto mitr
                = std::find(std::begin(mutations), std::end(mutations), new_mutation);
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
            auto extinct_mut = std::find(std::begin(mcounts), std::end(mcounts), 0);
            std::size_t mindex = mcounts.size();  // this will be the correct
                                                  // value if we use the else
                                                  // block below
            if (extinct_mut != std::end(mcounts)) // then we can recycle
                {
                    auto dist = std::distance(std::begin(mcounts), extinct_mut);
                    mutations[dist]
                        = std::move(new_mutation); // move new mutation into place
                    mindex = std::size_t(dist);    // update our value
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
        add_mutation_details(poptype &p, const std::vector<std::size_t> &mindexes,
                             const map_t &gams)
        {
            auto gam_recycling_bin = make_haploid_genome_queue(p.haploid_genomes);
            // Function object for calls to upper bound
            auto inserter
                = [&p](const double &__value, const std::size_t __mut) noexcept {
                      return __value < p.mutations[__mut].pos;
                  };

            for (const auto &gi : gams)
                {
                    auto n = p.haploid_genomes[gi.first].mutations;
                    auto s = p.haploid_genomes[gi.first].smutations;
                    for (auto mindex : mindexes)
                        {
                            auto pos = p.mutations[mindex].pos;

                            if (p.mutations[mindex].neutral)
                                {
                                    debug::validate_mutation_key_ranges(
                                        p.mutations, n.begin(), n.end());
                                    n.insert(std::upper_bound(n.begin(), n.end(), pos,
                                                              inserter),
                                             mindex);
                                }
                            else
                                {
                                    debug::validate_mutation_key_ranges(
                                        p.mutations, s.begin(), s.end());
                                    s.insert(std::upper_bound(s.begin(), s.end(), pos,
                                                              inserter),
                                             mindex);
                                }
                            // update mutation count
                            p.mcounts[mindex] +=
                                typename poptype::mcount_t::value_type(gi.second.size());
                        }
                    // get new haploid_genome
                    auto new_haploid_genome_key = recycle_haploid_genome(
                        p.haploid_genomes, gam_recycling_bin, n, s);
                    // update haploid_genome count
                    p.haploid_genomes[gi.first].n
                        -= decltype(p.haploid_genomes[gi.first].n)(gi.second.size());
                    p.haploid_genomes[new_haploid_genome_key].n += decltype(
                        p.haploid_genomes[new_haploid_genome_key].n)(gi.second.size());

                    // This updates every diploid to have a key to
                    // p.haploid_genomes[new_haploid_genome_key]
                    for (auto i : gi.second)
                        {
                            *i = new_haploid_genome_key;
                        }
                }
        }

        template <typename poptype,
                  typename
                  = std::enable_if<std::is_same<typename poptype::popmodel_t,
                                                fwdpp::poptypes::DIPLOID_TAG>::value>>
        std::unordered_map<std::size_t,
                           std::vector<typename poptype::diploid_type::first_type *>>
        collect_haploid_genomes(poptype &p, const std::vector<std::size_t> &indlist,
                                const std::vector<short> &clist)
        /*!
          \brief Helper function for add_mutation and add_mutations.

          Collects together all the haploid_genomes that will be updated, to minimize
          the number of
          new haploid_genomes created.
        */
        {
            std::unordered_map<std::size_t,
                               std::vector<typename poptype::diploid_type::first_type *>>
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
    } // namespace sugar

    template <typename poptype, class... Args>
    std::size_t
    add_mutation(poptype &p, const std::vector<std::size_t> &indlist,
                 const std::vector<short> &clist, Args &&... args)
    /*!
      \brief Add a mutation into a population at a given frequency.

      \param p A single locus population object.
      \param indlist A list of indexes of diploids into which to add the new
      mutations.
      \param clist A list of haploid_genomes.  See below.
      \param args Values required to cosnstruct a new mutation.  See below.

      \return The key referring to the location of the new mutation in the
      population

      Some notes:

      clist.size() must equal indlist.size()

      Values in \a clist must be 0, 1, or 2. These values mean to add the
      mutation to the first,
      second, or both haploid_genomes, resepectively, of each diploid in \a indlist.

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
                                   fwdpp::poptypes::DIPLOID_TAG>::value,
                      "poptype must be a single-locus object type");

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
                    throw std::out_of_range("clist contains elements < 0 and/or > 2");
            }
        if (indlist.size() != clist.size())
            throw std::invalid_argument("indlist and clist must be same length");

        // create a new mutation
        typename poptype::mutation_container::value_type new_mutant(
            std::forward<Args>(args)...);
        auto mindex = sugar::get_mut_index(p.mutations, p.mcounts, new_mutant);
        auto gams = sugar::collect_haploid_genomes(p, indlist, clist);
        sugar::add_mutation_details(p, {mindex}, gams);
        return mindex;
    }

    template <typename poptype>
    void
    add_mutations(poptype &p, const std::vector<std::size_t> &indlist,
                  const std::vector<short> &clist,
                  const std::vector<std::size_t> &mutation_indexes)
    /*!
      \brief Add a set of mutations into a set of individuals in a population.

	  \param p A diploid population object.
      \param indlist A list of indexes of diploids into which to add the new
      mutations.
      \param clist A list of haploid_genomes. See below.
      \param mutation_indexes The set of mutations to add.  See below.

      \return Nothing (void)

      Some notes:

      clist.size() must equal indlist.size()

      Values in \a clist must be 0, 1, or 2. These values mean to add the
      mutation to the first,
      second, or both haploid_genomes, resepectively, of each diploid in \a indlist.

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
                                   fwdpp::poptypes::DIPLOID_TAG>::value,
                      "poptype must be a diploid population ype");

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
                    throw std::invalid_argument(
                        "clist contains elements < 0 and/or > 2");
            }
        if (indlist.size() != clist.size())
            throw std::invalid_argument("indlist and clist must be same length");
        for (auto mi : mutation_indexes)
            {
                if (mi >= p.mutations.size())
                    throw std::out_of_range("mutation key out of range");
            }
        if (p.mcounts.size() != p.mutations.size())
            throw std::invalid_argument("p.mcounts.size() != p.mutations.size()");
        auto gams = sugar::collect_haploid_genomes(p, indlist, clist);
        sugar::add_mutation_details(p, mutation_indexes, gams);
    }

} // namespace fwdpp

#endif
