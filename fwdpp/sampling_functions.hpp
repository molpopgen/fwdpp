#ifndef __SAMPLING_FUNCTIONS_HPP__
#define __SAMPLING_FUNCTIONS_HPP__

#include <vector>
#include <map>
#include <functional>
#include <algorithm>
#include <utility>
#include <string>
#include <cassert>
#include <type_traits>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fwdpp/type_traits.hpp>
/*! @defgroup samplingPops Functions related to taking samples from simulated
  populations

  This collection of functions allows a user to draw a sample of size \f$n \ll
  2N\f$ from
  a simulated population.

  The library provides several overloads of the functions KTfwd::ms_sample and
  KTfwd::ms_sample_separate.

  The following features are in common to all versions of these functions:

  1. They only supp]ort bi-allelic mutation positions.
  2. All functions return some number of objects of type std::vector<
  std::pair<double, std::string> >.
  Each element in the vector corresponds to a mutation.  The double corresponds
  to the mutation position,
  and the string corresponds to the character states in each of the \f$n\f$
  samples.  The strings contain the character
  '0' when a sample has the ancestral state or '1' for the derived (mutant)
  state.  This type may be used
  to populate a polymorphism table from the
  [libsequence](http://molpopgen.github.io/libsequence/) library:
  \code
  auto x = KTfwd::ms_sample( appropriate arguments );
  Sequence::SimData xx(x.begin(),x.end());
  \endcode
  When a SimData object is written to a stream, its output format is the same
  as that used by Dick Hudson's coalescent simulation
  program [ms](http://home.uchicago.edu/~rhudson1/source/mksamples.html)
  3. All versions of KTfwd::ms_sample return vectors where the mutations
  affecting fitness and those not affecting fitness are intermingled.
  4. All versions of KTfwd::ms_sample_separate return pairs of vectors
  separating the mutations not affecting fitness from those that do.
  The first member of each pair is a vector of "neutral" mutations, and the
  second member is the vector of "selected" mutations.
  5. The functions are all implemented using sampling with replacement from the
  population.  Thus, setting \f$n = N\f$ (or \f$2N\f$) will NOT
  return all the variable sites in the entire population!!!

  A vector may contain data looking like the following:
  \verbatim
  0.0760 "10000"
  0.6805 "00010"
  \endverbatim

  The two rows are the two different positions (0.076 and 0.6805).  The first
  haplotype is "10", the second is "00", and the fourth is "01".

  If you converted the data to a Sequence::SimData object and printed it to
  screen or a file, the output would be:

  \verbatim
  //
  segsites: 2
  positions: 0.0760 0.6805
  10
  00
  00
  01
  00
  \endverbatim
*/

/*!
  @defgroup samplingPopsGamete Randomly-sampling individual gametes.
  @ingroup samplingPops

  These functions draws a sample of size \f$n \ll 2N\f$ from a simulated
  population.

  The gametes are randomly-sampled (with replacement) proportionally to their
  frequency in the population.

  The object passed to these functions is the container of gametes (e.g, some
  container of type KTfwd::gamete_base).

  The return values have no relation to any actual diploid individual in the
  population. Each haplotype
  corresponds to a simulated haplotype, but adjacent haplotypes are basically
  random draws of M&Ms from a jar.

  For these functions, \f$n\f$ can be odd or even.
*/

/*!
  @defgroup samplingPopsInd Randomly-sampling diploids.
  @ingroup samplingPops

  These functions pull a sample of \f$n\f$ simulated individuals from the
  population.

  The individuals are sampled uniformly, and with replacement, with no regard
  for their fitness.

  The object passed to these functions is a container of diploids.

  For these functions, \f$n\f$ must be an even number, as it represents the
  number of alleles to sample (twice the
  number of individuals).

  The return values of these functions store individuals in the order that they
  were sampled.
*/

namespace KTfwd
{
    /*!
      A variable site in a sample is a pair (pos,genotypes).

      This is equivalent to libsequence's Sequence::polymorphicSite
    */
    using sample_site_t = std::pair<double, std::string>;
    /*!
      A sample is a vector of variable sites.

      This is equivalent to libsequence's Sequence::polySiteVector.

      For this type, 'neutral' and 'selected' variants are intermingled.
    */
    using sample_t = std::vector<sample_site_t>;
    /*!
      A sample where 'neutral' and 'selected' variants are separated.

      'first' contains the 'neutral' variants, and 'second' contains the
      'selected' variants.
    */
    using sep_sample_t = std::pair<sample_t, sample_t>;

    /*!
      \brief Sampling from a population in an individual-based simulation
      \return A vector of variable sites
      \ingroup samplingPopsInd
    */
    template <typename mcont_t, typename gcont_t, typename allocator,
              typename diploid_geno_t,
              template <typename, typename> class vector_type>
    typename std::enable_if<traits::is_diploid<diploid_geno_t>::value,
                            sample_t>::type
    ms_sample(const gsl_rng *r, const mcont_t &mutations,
              const gcont_t &gametes,
              const vector_type<diploid_geno_t, allocator> &diploids,
              const unsigned &n, const bool &remove_fixed = true);

    /*!
      \brief Sampling from a population in an individual-based simulation.
      Selected and neutral mutations returned separately
      \return A pair of vectors of variable sites.  The first block is neutral
      variants, the second is non-neutral variants
      \ingroup samplingPopsInd
    */
    template <typename mcont_t, typename gcont_t, typename allocator,
              typename diploid_geno_t,
              template <typename, typename> class vector_type>
    typename std::enable_if<traits::is_diploid<diploid_geno_t>::value,
                            sep_sample_t>::type
    ms_sample_separate(const gsl_rng *r, const mcont_t &mutations,
                       const gcont_t &gametes,
                       const vector_type<diploid_geno_t, allocator> &diploids,
                       const unsigned &n, const bool &remove_fixed = true);

    /*!
      \brief Sample from an individual-based, multi-locus simulation.
      \return A vector of vectors of variable sites.  There is 1 vector per
      locus.
      \note Neutral + selected mutations intermixed
      \ingroup samplingPopsInd
    */
    template <typename mcont_t, typename gcont_t, typename dcont_t>
    typename std::enable_if<traits::is_diploid<typename dcont_t::value_type::
                                                   value_type>::value,
                            std::vector<sample_t>>::type
    ms_sample(const gsl_rng *r, const mcont_t &mutations,
              const gcont_t &gametes, const dcont_t &diploids,
              const unsigned &n, const bool &remove_fixed);

    /*!
      \brief Sample from an individual-based, multi-locus simulation.
      \return A vector of pairs of vectors of variable sites.  There is 1
      vector per locus.
      \note For each locus, the first member of the pair corresponds to neutral
      sites, the second to selected.
      \ingroup samplingPopsInd
    */
    template <typename mcont_t, typename gcont_t, typename dcont_t>
    typename std::enable_if<traits::is_diploid<typename dcont_t::value_type::
                                                   value_type>::value,
                            std::vector<sep_sample_t>>::type
    ms_sample_separate(const gsl_rng *r, const mcont_t &mutations,
                       const gcont_t &gametes, const dcont_t &diploids,
                       const unsigned &n, const bool &remove_fixed);
}
#endif
#include <fwdpp/sampling_functions.tcc>
