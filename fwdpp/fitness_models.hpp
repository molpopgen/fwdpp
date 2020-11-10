#ifndef _FITNESS_MODELS_HPP_
#define _FITNESS_MODELS_HPP_

#include <fwdpp/forward_types.hpp>
#include <fwdpp/fwd_functional.hpp>
#include <fwdpp/type_traits.hpp>
#include <fwdpp/util/named_type.hpp>
#include <cmath>
#include <stdexcept>
#include <type_traits>
#include <algorithm>
#include <functional>

/*
  Revision history

  KRT June 9 2016:
  1. Added range-based operator() to additive/multiplicative models.
  1a. I have not (yet) made these new operator() called by the old ones taking
  haploid_genomes as args.
  2. Updated documentation for site_dependent_genetic_value
  3. Changed typedef for result_type in additive/multiplicative models to refer
  to value from site_dependent_genetic_value

  KRT June 8 2016:
  1. Added range-based operator() to site_dependent_genetic_value.
  2. Changed site_dependent_genetic_value::operator() that took haploid_genomes as
  arguments
  to call the new range-based function.
  3. Fixed site_dependent_genetic_value::operator() for custom diploids.
  This had never been updated from the iterator-based library to the
  index-based library,
  but it was never caught b/c most calls go through the specific
  multiplicative/additive models.

  KRT April 24 2017
  1. additive_diploid and multiplicative_diploid changed to
  accept closures allowing them to calculate "fitness" or "trait value"
  more naturally.  This is in response to #49 on github.  This change
  does not break API compatibility.  The default is still to calculate
  "fitness".
*/

/*!
  \defgroup fitness Policies for calculating fitnesses and trait values.
  This group contains the following data structures to help you implement
  custom fitness policies:
  1. fwdpp::site_dependent_genetic_value for implementing typical
  population-genetic-like models where trait values are a function of the
  properties of individual mutations
  2. fwdpp::haplotype_dependent_fitness for implementing models where trait
  values are the functions properties of haplotypes.

  Truthfully, latter is so trivial that a library user may never see any need
  for it.

  The library also defines two site-dependent fitness models, which may also
  be used for trait value calculations:
  1. fwdpp::multiplicative_diploid
  2. fwdpp::additive_diploid

  These are arguably the "standard" models of the field as far as selection is
  concerned.

  Finally, fwdpp::no_selection is provided to force all diploid fitnesses to be
  equal to 1.
*/
namespace fwdpp
{
    /*! \brief Returns a fitness of 1
      \return A fitness of 1
      \param g1 A haploid_genome
      \param g2 A haploid_genome

      \note g1 and g2 must be part of the haploid_genome_base hierarchy
      \ingroup fitness
    */
    struct no_selection
    {
        /*!
          \brief Method for standard diploid simulations of a single locus.
        */
        using result_type = double;
        template <typename HaploidGenomeType, typename MutationContainerType>
        inline result_type
        operator()(const HaploidGenomeType &, const HaploidGenomeType &,
                   const MutationContainerType &) const noexcept
        {
            static_assert(traits::is_haploid_genome<HaploidGenomeType>::value,
                          "HaploidGenomeType::value_type must be a haploid_genome "
                          "type");
            static_assert(
                traits::is_mutation<typename MutationContainerType::value_type>::value,
                "MutationContainerType::value_type must be a mutation type");
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

    ///  \brief Function object for fitness/trait value as a function of
    ///  individual
    ///  mutations in a diploid
    ///
    ///  \note The updating policies must take a non-const reference to a
    ///  double
    ///  as the first argument and
    ///  an MutationContainerType::value_type as the second.  Further, they must not return
    ///  anything. Any remaining arguments needed should be passed via a
    ///  mechanism
    ///  such as std::bind and a function object, or via a lambda expression.
    ///  See
    ///  fwdpp::multiplicative_diploid for an example implementation.
    ///  \ingroup fitness
    struct site_dependent_genetic_value
    {
        //! The return value type
        using result_type = double;

        template <typename iterator_t, typename MutationContainerType,
                  typename updating_policy_hom, typename updating_policy_het,
                  typename make_return_value>
        inline result_type
        operator()(iterator_t first1, iterator_t last1, iterator_t first2,
                   iterator_t last2, const MutationContainerType &mutations,
                   const updating_policy_hom &fpol_hom,
                   const updating_policy_het &fpol_het,
                   const make_return_value &rv_function,
                   const double starting_value) const noexcept
        /*!
          Range-based call operator.  Calculates genetic values over ranges of
          mutation keys first1/last1
          and first2/last2, which are iterators derived from the 'smutations'
          of two haploid_genomes in a diploid.

          \param first1 Iterator to first mutation derived from haploid_genome 1
          \param last1 Iterator to one past the last mutation derived from
          haploid_genome 1
          \param first2 Iterator to first mutation derived from haploid_genome 2
          \param last2 Iterator to one past the last mutation derived from
          haploid_genome 2
          \param mutations The container of mutations for the simulation
          \param fpol_hom Policy that updates genetic value for a homozygous
          mutation.
          \param fpol_het Policy that updates genetic value for a heterozygous
          mutation.
          \param make_return_value Policy generated the final return value.
                 Must be equivalent to std::function<double(double)>
          \param starting_value The initial genetic value.

          \returns Fitness (double)
        */
        {
            static_assert(
                traits::is_mutation<typename MutationContainerType::value_type>::value,
                "MutationContainerType::value_type must be a mutation type");
            static_assert(
                std::is_convertible<
                    updating_policy_hom,
                    std::function<void(double &, const typename MutationContainerType::
                                                     value_type &)>>::value,
                "decltype(fpol_hom) must be convertible to "
                "std::function<void(double &,const typename "
                "MutationContainerType::value_type");
            static_assert(
                std::is_convertible<
                    updating_policy_het,
                    std::function<void(double &, const typename MutationContainerType::
                                                     value_type &)>>::value,
                "decltype(fpol_het) must be convertible to "
                "std::function<void(double &,const typename "
                "MutationContainerType::value_type");
            result_type w = starting_value;
            if (first1 == last1 && first2 == last2)
                return rv_function(w);
            else if (first1 == last1)
                {
                    for (; first2 != last2; ++first2)
                        fpol_het(w, mutations[*first2]);
                    return rv_function(w);
                }
            else if (first2 == last2)
                {
                    for (; first1 != last1; ++first1)
                        fpol_het(w, mutations[*first1]);
                    return rv_function(w);
                }
            for (; first1 != last1; ++first1)
                {
                    for (; first2 != last2 && *first1 != *first2
                           && mutations[*first2].pos < mutations[*first1].pos;
                         ++first2)
                        // All mutations in this range are Aa
                        {
                            fpol_het(w, mutations[*first2]);
                        }
                    if (first2 < last2
                        && (*first1 == *first2
                            || mutations[*first1].pos == mutations[*first2].pos))
                        // mutation with index first1 is homozygous
                        {
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
            return rv_function(w);
        }

        template <typename iterator_t, typename MutationContainerType,
                  typename updating_policy_hom, typename updating_policy_het>
        inline result_type
        operator()(iterator_t first1, iterator_t last1, iterator_t first2,
                   iterator_t last2, const MutationContainerType &mutations,
                   const updating_policy_hom &fpol_hom,
                   const updating_policy_het &fpol_het,
                   const double starting_value) const noexcept
        {
            return this->operator()(
                first1, last1, first2, last2, mutations, fpol_hom, fpol_het,
                [](double d) { return d; }, starting_value);
        }

        template <typename HaploidGenomeType, typename MutationContainerType,
                  typename updating_policy_hom, typename updating_policy_het,
                  typename make_return_value>
        inline result_type
        operator()(const HaploidGenomeType &g1, const HaploidGenomeType &g2,
                   const MutationContainerType &mutations,
                   const updating_policy_hom &fpol_hom,
                   const updating_policy_het &fpol_het,
                   const make_return_value &rv_function,
                   const double starting_value) const noexcept
        /*!
          Calculates genetic value for a diploid whose genotype
          across sites is given by haploid_genomes g1 and  g2.

          \param g1 A haploid_genome
          \param g2 A haploid_genome
          \param mutations The container of mutations for the simulation
          \param fpol_hom Policy that updates genetic value for a homozygous
          mutation.
          \param fpol_het Policy that updates genetic value for a heterozygous
          mutation.
          \param starting_value The initial genetic value.

          \returns Fitness (double)
        */
        {
            static_assert(traits::is_haploid_genome<HaploidGenomeType>::value,
                          "HaploidGenomeType::value_type must be a haploid_genome "
                          "type");
            return this->operator()(g1.smutations.cbegin(), g1.smutations.cend(),
                                    g2.smutations.cbegin(), g2.smutations.cend(),
                                    mutations, fpol_hom, fpol_het, rv_function,
                                    starting_value);
        }

        template <typename HaploidGenomeType, typename MutationContainerType,
                  typename updating_policy_hom, typename updating_policy_het>
        inline result_type
        operator()(const HaploidGenomeType &g1, const HaploidGenomeType &g2,
                   const MutationContainerType &mutations,
                   const updating_policy_hom &fpol_hom,
                   const updating_policy_het &fpol_het,
                   const double starting_value) const noexcept
        {
            return this->operator()(
                g1, g2, mutations, fpol_hom, fpol_het, [](double d) { return d; },
                starting_value);
        }

        template <typename DiploidType, typename GenomeContainerType,
                  typename MutationContainerType, typename updating_policy_hom,
                  typename updating_policy_het, typename make_return_value>
        inline result_type
        operator()(const DiploidType &dip, const GenomeContainerType &haploid_genomes,
                   const MutationContainerType &mutations,
                   const updating_policy_hom &fpol_hom,
                   const updating_policy_het &fpol_het,
                   const make_return_value &rv_function,
                   const double starting_value) const noexcept
        /*!
          Calculates genetic value for a diploid type.

          \param dip A diploid
          \param haploid_genomes The container of haploid_genomes for the simulation.
          \param mutations The container of mutations for the simulation
          \param fpol_hom Policy that updates genetic value for a homozygous
          mutation.
          \param fpol_het Policy that updates genetic value for a heterozygous
          mutation.
          \param starting_value The initial genetic value.

          \returns Fitness (double)
        */
        {
            return this->operator()(haploid_genomes[dip.first],
                                    haploid_genomes[dip.second], mutations, fpol_hom,
                                    fpol_het, rv_function, starting_value);
        }

        template <typename DiploidType, typename GenomeContainerType,
                  typename MutationContainerType, typename updating_policy_hom,
                  typename updating_policy_het>
        inline result_type
        operator()(const DiploidType &dip, const GenomeContainerType &haploid_genomes,
                   const MutationContainerType &mutations,
                   const updating_policy_hom &fpol_hom,
                   const updating_policy_het &fpol_het,
                   const double starting_value) const noexcept
        {
            return this->operator()(
                dip, haploid_genomes, mutations, fpol_hom, fpol_het,
                [](double d) { return d; }, starting_value);
        }
    };

    /// Typedef for API compatibility.
    using site_dependent_fitness = site_dependent_genetic_value;

    /// \brief Function object for fitness/trait value
    ///  as a function of the two haplotypes
    ///  in a diploid
    ///
    /// \note This really is just a convenience function. Depending on the
    /// specifics of the model, this function may be totally unnecessary.
    /// \ingroup fitness
    struct haplotype_dependent_trait_value
    {
        using result_type = double;
        template <typename HaploidGenomeType, typename MutationContainerType,
                  typename haplotype_policy, typename DiploidTypePolicy>
        inline result_type
        operator()(const HaploidGenomeType &g1, const HaploidGenomeType &g2,
                   const MutationContainerType &mutations, const haplotype_policy &hpol,
                   const DiploidTypePolicy &dpol) const noexcept
        ///
        ///\param g1 A haploid_genome
        ///\param g2 A haploid_genome
        ///\param mutations Container of mutations
        ///\param hpol A policy whose first argument is an iterator to
        /// a haploid_genome.
        /// Remaining arguments may be bound via std::bind or the
        /// equivalent.  The
        /// policy returns a double representing the effect of this
        /// haplotype on
        /// genetic value
        ///\param dpol A policy whose first two arguments are doubles,
        /// each of which
        /// represents the effect of g1 and g2, respectively.  Remaining
        /// arguments
        /// may be bound via std::bind or the equivalent.  The policy
        /// returns a
        /// double representing the genetic value of a diploid genotype g1/g2
        ///\return dpol( hpol(g1,mutations), hpol(g2,mutations) )
        ///
        {
            static_assert(traits::is_haploid_genome_v<HaploidGenomeType>,
                          "HaploidGenomeType must be a haploid_genome type");
            static_assert(
                traits::is_mutation<typename MutationContainerType::value_type>::value,
                "MutationContainerType::value_type must be a mutation type");
            return dpol(hpol(g1, mutations), hpol(g2, mutations));
        }
        template <typename DiploidType, typename GenomeContainerType,
                  typename MutationContainerType, typename haplotype_policy,
                  typename DiploidTypePolicy>
        inline result_type
        operator()(const DiploidType &diploid,
                   const GenomeContainerType &haploid_genomes,
                   const MutationContainerType &mutations, const haplotype_policy &hpol,
                   const DiploidTypePolicy &dpol) const noexcept
        ///\param diploid a diploid
        ///\param haploid_genomes Container of haploid_genomes
        ///\param mutations Container of mutations
        ///\param hpol A policy whose first argument is an iterator to a
        /// haploid_genome.
        /// Remaining arguments may be bound via std::bind or the
        /// equivalent.  The
        /// policy returns a double representing the effect of this
        /// haplotype on
        /// genetic value
        ///\param dpol A policy whose first two arguments are doubles,
        /// each of which
        /// represents the effect of g1 and g2, respectively.  Remaining
        /// arguments
        /// may be bound via std::bind or the equivalent.  The policy
        /// returns a
        /// double representing the genetic value of a diploid genotype g1/g2
        ///\return dpol( hpol(g1,mutations), hpol(g2,mutations) )
        {
            static_assert(traits::is_diploid<DiploidType>::value,
                          "diploid_t must represent a diploid");
            return this->operator()(haploid_genomes[diploid.first],
                                    haploid_genomes[diploid.second], mutations, hpol,
                                    dpol);
        }
    };

    /// Typedef for backwards API compatibility.
    using haplotype_dependent_fitness = haplotype_dependent_trait_value;

    struct genetic_value_is_trait
    {
    };

    struct genetic_value_is_fitness
    {
    };

    /// Strong wrapper around a double signifying the scaling of a model.
    /// Typical use would be to assign values 0, sh, scaling*s to
    /// genotypes AA, Aa, aa
    using trait = strong_types::named_type<double, genetic_value_is_trait>;

    /// Strong wrapper around a double signifying the scaling of a model.
    /// Typical use would be to assign values 0, sh, scaling*s to
    /// genotypes AA, Aa, aa
    using fitness = strong_types::named_type<double, genetic_value_is_fitness>;

    inline bool
    assign_is_trait_value(const trait &)
    /// Helper function for genetic value object constructors
    {
        return true;
    }

    inline bool
    assign_is_trait_value(const fitness &)
    /// Helper function for genetic value object constructors
    {
        return false;
    }

    inline bool
    assign_is_fitness_value(const trait &)
    /// Helper function for genetic value object constructors
    {
        return false;
    }

    inline bool
    assign_is_fitness_value(const fitness &)
    /// Helper function for genetic value object constructors
    {
        return true;
    }

    /// \brief Multiplicative fitness or trait value across sites
    /// This function object calculate the genetic
    /// value of a diploid according to an multiplicative model
    /// with dominance effects per mutation.
    ///
    /// The genetic value is calculated via
    /// fwdpp::site_dependent_genetic_value.
    ///
    /// The genetic value may be treated as a fitness
    /// or as a trait value based on arguments passed
    /// to the constructor.
    ///
    /// \ingroup fitness
    struct multiplicative_diploid
    {
        std::function<double(double)>
        assign_f(trait &)
        {
            return [](const double d) { return d - 1.0; };
        }
        std::function<double(double)>
        assign_f(fitness &)
        {
            return [](const double d) { return std::max(0.0, d); };
        }
        const double scaling;
        const bool gvalue_is_trait;
        const bool gvalue_is_fitness;
        using result_type = site_dependent_genetic_value::result_type;
        const std::function<double(double)> make_return_value;
        multiplicative_diploid(trait t)
            : scaling{t.get()}, gvalue_is_trait(assign_is_trait_value(t)),
              gvalue_is_fitness(assign_is_fitness_value(t)), make_return_value{
                                                                 assign_f(t)}
        /// Construct an object to calculate trait/phenotype values.
        /// \param t fwdpp::trait, where the double repsresents the scaling of "aa" trait values.
        {
            if (!std::isfinite(scaling))
                {
                    throw std::invalid_argument("scaling parameter must be finite");
                }
        }

        template <typename make_return_value_fxn>
        multiplicative_diploid(fitness f, make_return_value_fxn &&make_rv)
            : scaling{f.get()}, gvalue_is_trait{assign_is_trait_value(f)},
              gvalue_is_fitness{assign_is_fitness_value(f)},
              make_return_value{std::forward<make_return_value_fxn>(make_rv)}
        {
            if (!std::isfinite(scaling))
                {
                    throw std::invalid_argument("scaling parameter must be finite");
                }
        }

        multiplicative_diploid(fitness gvtype)
            : scaling{gvtype.get()}, gvalue_is_trait(assign_is_trait_value(gvtype)),
              gvalue_is_fitness(assign_is_fitness_value(gvtype)), make_return_value{
                                                                      assign_f(gvtype)}
        /// Construct an object to calculate fitness values.
        /// \param gvtype fwdpp::fitness, where the double repsresents the scaling of "aa" trait values.
        {
            if (!std::isfinite(scaling))
                {
                    throw std::invalid_argument("scaling parameter must be finite");
                }
        }

        template <typename make_return_value_fxn>
        multiplicative_diploid(trait t, make_return_value_fxn &&make_rv)
            : scaling{t.get()}, gvalue_is_trait{assign_is_trait_value(t)},
              gvalue_is_fitness{assign_is_fitness_value(t)},
              make_return_value{std::forward<make_return_value_fxn>(make_rv)}
        {
            if (!std::isfinite(scaling))
                {
                    throw std::invalid_argument("scaling parameter must be finite");
                }
        }

        template <typename iterator_t, typename MutationContainerType>
        inline result_type
        operator()(iterator_t first1, iterator_t last1, iterator_t first2,
                   iterator_t last2, const MutationContainerType &mutations) const
            noexcept
        /// Range-based overload
        {
            using __mtype = typename MutationContainerType::value_type;
            return site_dependent_genetic_value()(
                first1, last1, first2, last2, mutations,
                [this](double &value, const __mtype &mut) noexcept {
                    value *= (1. + scaling * mut.s);
                },
                [](double &value, const __mtype &mut) noexcept {
                    value *= (1. + mut.h * mut.s);
                },
                make_return_value, 1.);
        }

        template <typename HaploidGenomeType, typename MutationContainerType>
        inline double
        operator()(const HaploidGenomeType &g1, const HaploidGenomeType &g2,
                   const MutationContainerType &mutations) const noexcept
        ///  \param g1 A haploid_genome
        ///  \param g2 A haploid_genome
        ///  \param mutations Container of mutations
        ///  \return Multiplicative genetic value across sites.
        {
            using __mtype = typename MutationContainerType::value_type;
            return site_dependent_genetic_value()(
                g1, g2, mutations,
                [this](double &value, const __mtype &mut) noexcept {
                    value *= (1. + scaling * mut.s);
                },
                [](double &value, const __mtype &mut) noexcept {
                    value *= (1. + mut.h * mut.s);
                },
                make_return_value, 1.);
        }

        template <typename DiploidType, typename GenomeContainerType,
                  typename MutationContainerType>
        inline result_type
        operator()(const DiploidType &dip, const GenomeContainerType &haploid_genomes,
                   const MutationContainerType &mutations) const noexcept
        ///  \brief Overload for diploids.  This is what a programmer's
        ///  functions will call.
        {
            return this->operator()(haploid_genomes[dip.first],
                                    haploid_genomes[dip.second], mutations);
        }
    };

    /// \brief Additive fitness or trait value across sites
    /// This function object calculate the genetic
    /// value of a diploid according to an additive model
    /// with dominance effects per mutation.
    ///
    /// The genetic value is calculated via
    /// fwdpp::site_dependent_genetic_value.
    ///
    /// The genetic value may be treated as a fitness
    /// or as a trait value based on arguments passed
    /// to the constructor.
    ///
    /// \ingroup fitness
    struct additive_diploid
    {
        /// Specifies final mapping of genetic value
        std::function<double(double)> assign_f(trait)
        {
            return [](const double d) { return d; };
        }
        std::function<double(double)> assign_f(fitness)
        {
            return [](const double d) { return std::max(0.0, 1.0 + d); };
        }
        const double scaling;
        const bool gvalue_is_trait;
        const bool gvalue_is_fitness;
        const std::function<double(double)> make_return_value;
        using result_type = site_dependent_genetic_value::result_type;

        additive_diploid(fitness f)
            : scaling{f.get()}, gvalue_is_trait{assign_is_trait_value(f)},
              gvalue_is_fitness{assign_is_fitness_value(f)}, make_return_value{
                                                                 assign_f(f)}
        /// Construct an object to calculate fitness values
        /// \param f fwdpp::fitness, where the double represents the scaling of the "aa" genotype
        {
            if (!std::isfinite(scaling))
                {
                    throw std::invalid_argument("scaling parameter must be finite");
                }
        }

        template <typename make_return_value_fxn>
        additive_diploid(fitness f, make_return_value_fxn &&make_rv)
            : scaling{f.get()}, gvalue_is_trait{assign_is_trait_value(f)},
              gvalue_is_fitness{assign_is_fitness_value(f)},
              make_return_value{std::forward<make_return_value_fxn>(make_rv)}
        {
            if (!std::isfinite(scaling))
                {
                    throw std::invalid_argument("scaling parameter must be finite");
                }
        }

        additive_diploid(trait t)
            : scaling{t.get()}, gvalue_is_trait{assign_is_trait_value(t)},
              gvalue_is_fitness{assign_is_fitness_value(t)}, make_return_value{
                                                                 assign_f(t)}
        /// Construct an object to calculate trait/phenotype values
        /// \param t fwdpp::trait, where the double represents the scaling of the "aa" genotype
        {
            if (!std::isfinite(scaling))
                {
                    throw std::invalid_argument("scaling parameter must be finite");
                }
        }

        template <typename make_return_value_fxn>
        additive_diploid(trait t, make_return_value_fxn &&make_rv)
            : scaling{t.get()}, gvalue_is_trait{assign_is_trait_value(t)},
              gvalue_is_fitness{assign_is_fitness_value(t)},
              make_return_value{std::forward<make_return_value_fxn>(make_rv)}
        {
            if (!std::isfinite(scaling))
                {
                    throw std::invalid_argument("scaling parameter must be finite");
                }
        }

        template <typename iterator_t, typename MutationContainerType>
        inline result_type
        operator()(iterator_t first1, iterator_t last1, iterator_t first2,
                   iterator_t last2, const MutationContainerType &mutations) const
            noexcept
        /// Range-based overload
        {
            using __mtype = typename MutationContainerType::value_type;
            return site_dependent_genetic_value()(
                first1, last1, first2, last2, mutations,
                [this](double &value, const __mtype &mut) noexcept {
                    value += (scaling * mut.s);
                },
                [](double &value, const __mtype &mut) noexcept {
                    value += (mut.h * mut.s);
                },
                make_return_value, 0.);
        }

        template <typename HaploidGenomeType, typename MutationContainerType>
        inline result_type
        operator()(const HaploidGenomeType &g1, const HaploidGenomeType &g2,
                   const MutationContainerType &mutations) const noexcept
        ///  \param g1 A haploid_genome
        ///  \param g2 A haploid_genome
        ///  \param mutations A container of mutations
        ///  \return Additive genetic value across sites.
        ///  \note g1 and g2 must be part of the haploid_genome_base hierarchy
        {
            using __mtype = typename MutationContainerType::value_type;
            return site_dependent_genetic_value()(
                g1, g2, mutations,
                [this](double &value, const __mtype &mut) noexcept {
                    value += (scaling * mut.s);
                },
                [](double &value, const __mtype &mut) noexcept {
                    value += (mut.h * mut.s);
                },
                make_return_value, 0.);
        }

        ///  \brief Overload for diploids.  This is what a programmer's
        ///  functions will call.
        template <typename DiploidType, typename GenomeContainerType,
                  typename MutationContainerType>
        inline result_type
        operator()(const DiploidType &dip, const GenomeContainerType &haploid_genomes,
                   const MutationContainerType &mutations) const noexcept
        {
            return this->operator()(haploid_genomes[dip.first],
                                    haploid_genomes[dip.second], mutations);
        }
    };
} // namespace fwdpp
#endif /* _FITNESS_MODELS_HPP_ */
