#ifndef _FITNESS_MODELS_HPP_
#define _FITNESS_MODELS_HPP_

#include <fwdpp/forward_types.hpp>
#include <fwdpp/fwd_functional.hpp>
#include <fwdpp/type_traits.hpp>
#include <cassert>
#include <type_traits>
#include <algorithm>
#include <functional>

/*
  Revision history

  KRT June 9 2016:
  1. Added range-based operator() to additive/multiplicative models.
  1a. I have not (yet) made these new operator() called by the old ones taking
  gametes as args.
  2. Updated documentation for site_dependent_fitness
  3. Changed typedef for result_type in additive/multiplicative models to refer
  to value from site_dependent fitness

  KRT June 8 2016:
  1. Added range-based operator() to site_dependent_fitness.
  2. Changed site_dependent_fitness::operator() that took gametes as arguments
  to call the new range-based function.
  3. Fixed site_dependent_fitness::operator() for custom diploids.
  This had never been updated from the iterator-based library to the
  index-based library,
  but it was never caught b/c most calls go through the specific
  multiplicative/additive models.
*/

/*!
  \defgroup fitness Policies for calculating fitnesses.
  This group contains the following data structures to help you implement
  custom fitness policies:
  1. KTfwd::site_dependent_fitness for implementing typical
  population-genetic-like models where trait values are a function of the
  properties of individual mutations
  2. KTfwd::haplotype_dependent_fitness for implementing models where trait
  values are the functions properties of haplotypes.

  Truthfully, latter is so trivial that a library user may never see any need
  for it.

  The library also defines two site-dependent fitness models:
  1. KTfwd::multiplicative_diploid
  2. KTfwd::additive_diploid

  These are arguably the "standard" models of the field as far as selection is
  concerned.

  Finally, KTfwd::no_selection is provided to force all diploid fitnesses to be
  equal to 1.
*/
namespace KTfwd
{
    /*! \brief Returns a fitness of 1
      \return A fitness of 1
      \param g1 A gamete
      \param g2 A gamete

      \note g1 and g2 must be part of the gamete_base hierarchy
      \ingroup fitness
    */
    struct no_selection
    {
        /*!
          \brief Method for standard diploid simulations of a single locus.
        */
        using result_type = double;
        template <typename gamete_type, typename mcont_t>
        inline result_type
        operator()(const gamete_type &, const gamete_type &,
                   const mcont_t &) const noexcept
        {
            static_assert(traits::is_gamete<gamete_type>::value,
                          "gamete_type::value_type must be a gamete type");
            static_assert(
                traits::is_mutation<typename mcont_t::value_type>::value,
                "mcont_t::value_type must be a mutation type");
            return 1.;
        }
        //! \brief Naive implementation for non-standard cases
        template <typename T>
        inline result_type
        operator()(const T &) const noexcept
        {
            return 1.;
        }
    };

    /*! \brief Function object for fitness as a function of individual
      mutations in a diploid
      \note The updating policies must take a non-const reference to a double
      as the first argument and
      an mcont_t::value_type as the second.  Further, they must not return
      anything. Any remaining arguments needed should be passed via a mechanism
      such as std::bind and a function object, or via a lambda expression.  See
      KTfwd::multiplicative_diploid for an example implementation.
      \ingroup fitness
    */
    struct site_dependent_fitness
    {
        //! The return value type
        using result_type = double;

        template <typename iterator_t, typename mcont_t,
                  typename fitness_updating_policy_hom,
                  typename fitness_updating_policy_het>
        inline result_type
        operator()(iterator_t first1, iterator_t last1, iterator_t first2,
                   iterator_t last2, const mcont_t &mutations,
                   const fitness_updating_policy_hom &fpol_hom,
                   const fitness_updating_policy_het &fpol_het,
                   const double &starting_fitness = 1.) const noexcept
        /*!
          Range-based call operator.  Calculates fitnesses over ranges of
          mutation keys first1/last1
          and first2/last2, which are iterators derived from the 'smutations'
          of two gametes in a diploid.

          \param first1 Iterator to first mutation derived from gamete 1
          \param last1 Iterator to one past the last mutation derived from
          gamete 1
          \param first2 Iterator to first mutation derived from gamete 2
          \param last2 Iterator to one past the last mutation derived from
          gamete 2
          \param mutations The container of mutations for the simulation
          \param fpol_hom Policy that updates fitness for a homozygous
          mutation.
          \param fpol_het Policy that updates fitness for a heterozygous
          mutation.
          \param starting_fitness The initial fitness value.

          \returns Fitness (double)
        */
        {
            static_assert(
                traits::is_mutation<typename mcont_t::value_type>::value,
                "mcont_t::value_type must be a mutation type");
            static_assert(
                std::is_convertible<fitness_updating_policy_hom,
                                    std::function<void(
                                        double &, const typename mcont_t::
                                                      value_type &)>>::value,
                "decltype(fpol_hom) must be convertible to "
                "std::function<void(double &,const typename "
                "mcont_t::value_type");
            static_assert(
                std::is_convertible<fitness_updating_policy_het,
                                    std::function<void(
                                        double &, const typename mcont_t::
                                                      value_type &)>>::value,
                "decltype(fpol_het) must be convertible to "
                "std::function<void(double &,const typename "
                "mcont_t::value_type");
            result_type w = starting_fitness;
            if (first1 == last1 && first2 == last2)
                return w;
            else if (first1 == last1)
                {
                    for (; first2 != last2; ++first2)
                        fpol_het(w, mutations[*first2]);
                    return w;
                }
            else if (first2 == last2)
                {
                    for (; first1 != last1; ++first1)
                        fpol_het(w, mutations[*first1]);
                    return w;
                }
            for (; first1 != last1; ++first1)
                {
                    for (;
                         first2 != last2 && *first1 != *first2
                         && !(mutations[*first2].pos > mutations[*first1].pos);
                         ++first2)
                        // All mutations in this range are Aa
                        {
                            assert(mutations[*first2].pos
                                   < mutations[*first1].pos);
                            fpol_het(w, mutations[*first2]);
                        }
                    if (first2 < last2 && *first1 == *first2) // mutation with
                                                              // index first1
                                                              // is homozygous
                        {
                            assert(mutations[*first2].pos
                                   == mutations[*first1].pos);
                            fpol_hom(w, mutations[*first1]);
                            ++first2; // increment so that we don't re-process
                                      // this site as a het next time 'round
                        }
                    else // mutation first1 is heterozygous
                        {
                            fpol_het(w, mutations[*first1]);
                        }
                }
            for (; first2 != last2; ++first2)
                {
                    fpol_het(w, mutations[*first2]);
                }
            return w;
        }

        template <typename gamete_type, typename mcont_t,
                  typename fitness_updating_policy_hom,
                  typename fitness_updating_policy_het>
        inline result_type
        operator()(const gamete_type &g1, const gamete_type &g2,
                   const mcont_t &mutations,
                   const fitness_updating_policy_hom &fpol_hom,
                   const fitness_updating_policy_het &fpol_het,
                   const double &starting_fitness = 1.) const noexcept
        /*!
          Calculates fitnesses for a diploid whose genotype
          across sites is given by gametes g1 and  g2.

          \param g1 A gamete
          \param g2 A gamete
          \param mutations The container of mutations for the simulation
          \param fpol_hom Policy that updates fitness for a homozygous
          mutation.
          \param fpol_het Policy that updates fitness for a heterozygous
          mutation.
          \param starting_fitness The initial fitness value.

          \returns Fitness (double)
        */
        {
            static_assert(traits::is_gamete<gamete_type>::value,
                          "gamete_type::value_type must be a gamete type");
            return this->operator()(
                g1.smutations.cbegin(), g1.smutations.cend(),
                g2.smutations.cbegin(), g2.smutations.cend(), mutations,
                fpol_hom, fpol_het, starting_fitness);
        }

        template <typename diploid2dispatch, typename gcont_t,
                  typename mcont_t, typename fitness_updating_policy_hom,
                  typename fitness_updating_policy_het>
        inline result_type
        operator()(const diploid2dispatch &dip, const gcont_t &gametes,
                   const mcont_t &mutations,
                   const fitness_updating_policy_hom &fpol_hom,
                   const fitness_updating_policy_het &fpol_het,
                   const double &starting_fitness = 1.) const noexcept
        /*!
          Calculates fitnesses for a custom diploid type. See @ref
          md_md_customdip.

          \param dip A diploid
          \param gametes The container of gametes for the simulation.
          \param mutations The container of mutations for the simulation
          \param fpol_hom Policy that updates fitness for a homozygous
          mutation.
          \param fpol_het Policy that updates fitness for a heterozygous
          mutation.
          \param starting_fitness The initial fitness value.

          \returns Fitness (double)
        */
        {
            return this->operator()(gametes[dip.first], gametes[dip.second],
                                    mutations, fpol_hom, fpol_het,
                                    starting_fitness);
        }
    };

    /*! \brief Function object for fitness as a function of the two haplotypes
      in a diploid

      \param g1 A gamete
      \param g2 A gamete
      \param mutation Container of mutations
      \param hpol A policy whose first argument is an iterator to a gamete.
      Remaining arguments may be bound via std::bind or the equivalent.  The
      policy returns a double representing the effect of this haplotype on
      fitness
      \param dpol A policy whose first two arguments are doubles, each of which
      represents the effect of g1 and g2, respectively.  Remaining arguments
      may be bound via std::bind or the equivalent.  The policy returns a
      double representing the fitness of a diploid genotype g1/g2
      \return dpol( hpol(g1,mutations), hpol(g2,mutations) )
      \note This really is just a convenience function. Depending on the
      specifics of the model, this function may be totally unnecessary.
      \ingroup fitness
    */
    struct haplotype_dependent_fitness
    {
        using result_type = double;
        template <typename gamete_type, typename mcont_t,
                  typename haplotype_policy, typename diploid_policy>
        inline result_type
        operator()(const gamete_type &g1, const gamete_type &g2,
                   const mcont_t &mutations, const haplotype_policy &hpol,
                   const diploid_policy &dpol) const noexcept
        {
            static_assert(typename traits::is_gamete<gamete_type>::type(),
                          "gamete_type must be a gamete type");
            static_assert(
                traits::is_mutation<typename mcont_t::value_type>::value,
                "mcont_t::value_type must be a mutation type");
            return dpol(hpol(g1, mutations), hpol(g2, mutations));
        }
        template <typename diploid_t, typename gcont_t, typename mcont_t,
                  typename haplotype_policy, typename diploid_policy>
        inline result_type
        operator()(const diploid_t &diploid, const gcont_t &gametes,
                   const mcont_t &mutations, const haplotype_policy &hpol,
                   const diploid_policy &dpol) const noexcept
        {
            static_assert(traits::is_diploid<diploid_t>::value,
                          "diploid_t must represent a diploid");
            return this->operator()(gametes[diploid.first],
                                    gametes[diploid.second], mutations, hpol,
                                    dpol);
        }
    };

    /*! \brief Multiplicative fitness across sites
      \param g1 A gamete
      \param g2 A gamete
      \param mutation Container of mutations
      \param scaling Fitnesses are 1, 1+h*s, 1+scaling*s, for AA,Aa,aa, resp.
      This parameter lets you make sure your
      simulation is on the same scale as various formula in the literature
      \return Multiplicative fitness across sites.
      \ingroup fitness
    */
    struct multiplicative_diploid
    {
        using result_type = site_dependent_fitness::result_type;
        template <typename iterator_t, typename mcont_t>
        inline result_type
        operator()(iterator_t first1, iterator_t last1, iterator_t first2,
                   iterator_t last2, const mcont_t &mutations,
                   const double &scaling = 1.) const noexcept
        {
            using __mtype = typename mcont_t::value_type;
            return std::max(
                0., site_dependent_fitness()(
                        first1, last1, first2, last2, mutations,
                        [&](double &fitness, const __mtype &mut) noexcept {
                            fitness *= (1. + scaling * mut.s);
                        },
                        [](double &fitness, const __mtype &mut) noexcept {
                            fitness *= (1. + mut.h * mut.s);
                        },
                        1.));
        }

        template <typename gamete_type, typename mcont_t>
        inline double
        operator()(const gamete_type &g1, const gamete_type &g2,
                   const mcont_t &mutations, const double &scaling = 1.) const
            noexcept
        {
            using __mtype = typename mcont_t::value_type;
            return std::max(
                0., site_dependent_fitness()(
                        g1, g2, mutations,
                        [&](double &fitness, const __mtype &mut) noexcept {
                            fitness *= (1. + scaling * mut.s);
                        },
                        [](double &fitness, const __mtype &mut) noexcept {
                            fitness *= (1. + mut.h * mut.s);
                        },
                        1.));
        }
        /*!
          \brief Overload for custom diploids.  This is what a programmer's
          functions will call.
          \note See @ref md_md_customdip
        */
        template <typename diploid2dispatch, typename gcont_t,
                  typename mcont_t>
        inline result_type
        operator()(const diploid2dispatch &dip, const gcont_t &gametes,
                   const mcont_t &mutations, const double &scaling = 1.) const
            noexcept
        {
            return this->operator()(gametes[dip.first], gametes[dip.second],
                                    mutations, scaling);
        }
    };

    /*! \brief Additive fitness across sites
      \param g1 A gamete
      \param g2 A gamete
      \param mutations A container of mutations
      \param scaling Fitnesses are 1, 1+h*s, 1+scaling*s, for AA,Aa,aa, resp.
      This parameter lets you make sure your
      simulation is on the same scale as various formula in the literature
      \return Additive fitness across sites.
      \note g1 and g2 must be part of the gamete_base hierarchy
      \ingroup fitness
    */
    struct additive_diploid
    {
        using result_type = site_dependent_fitness::result_type;
        template <typename iterator_t, typename mcont_t>
        inline result_type
        operator()(iterator_t first1, iterator_t last1, iterator_t first2,
                   iterator_t last2, const mcont_t &mutations,
                   const double &scaling = 1.) const noexcept
        {
            using __mtype = typename mcont_t::value_type;
            return std::max(
                0.,
                1. + site_dependent_fitness()(
                         first1, last1, first2, last2, mutations,
                         [&](double &fitness, const __mtype &mut) noexcept {
                             fitness += (scaling * mut.s);
                         },
                         [](double &fitness, const __mtype &mut) noexcept {
                             fitness += (mut.h * mut.s);
                         },
                         0.));
        }

        template <typename gamete_type, typename mcont_t>
        inline result_type
        operator()(const gamete_type &g1, const gamete_type &g2,
                   const mcont_t &mutations, const double &scaling = 1.) const
            noexcept
        {
            using __mtype = typename mcont_t::value_type;
            return std::max(
                0.,
                1. + site_dependent_fitness()(
                         g1, g2, mutations,
                         [=](double &fitness, const __mtype &mut) noexcept {
                             fitness += (scaling * mut.s);
                         },
                         [](double &fitness, const __mtype &mut) noexcept {
                             fitness += (mut.h * mut.s);
                         },
                         0.));
        }

        /*!
          \brief Overload for custom diploids.  This is what a programmer's
          functions will call.
          \note See @ref md_md_customdip
        */
        template <typename diploid2dispatch, typename gcont_t,
                  typename mcont_t>
        inline result_type
        operator()(const diploid2dispatch &dip, const gcont_t &gametes,
                   const mcont_t &mutations, const double &scaling = 1.) const
            noexcept
        {
            return this->operator()(gametes[dip.first], gametes[dip.second],
                                    mutations, scaling);
        }
    };
}
#endif /* _FITNESS_MODELS_HPP_ */
