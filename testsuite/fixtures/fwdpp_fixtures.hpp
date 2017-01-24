/*!
  \file fwdpp_fixtures.hpp
  \brief Minimal fixtures for testing main library with lowest-level (e.g.
  non-sugar) types
  \ingroup unit
*/
#ifndef FWDPP_TESTSUITE_FWDPP_FIXTURES_HPP
#define FWDPP_TESTSUITE_FWDPP_FIXTURES_HPP

#include <vector>
#include <utility>
#include <unordered_set>
#include <fwdpp/fwd_functional.hpp>
#include <fwdpp/forward_types.hpp>
#include <gsl/gsl_rng.h>

/*
 * These fixtures are for the low-level fwdpp
 * API.  Thus, we need a bunch of typedefs
 * to maintain our sanity. "Real-world" use
 * of fwdpp make use of "sugar" types (fwdpp/sugar/*.hpp),
 * which handles a lot of this stuff for you.
 */
using mtype = KTfwd::mutation;
using mcont_t = std::vector<mtype>;
using gcont_t = std::vector<KTfwd::gamete>;
using dipvector_t = std::vector<std::pair<std::size_t, std::size_t>>;
using lookup_table_t
    = std::unordered_set<double, std::hash<double>, KTfwd::equal_eps>;
using mcounts_t = std::vector<KTfwd::uint_t>;

struct standard_empty_single_deme_fixture
/*!
  Basic stuff needed for a simulation of a single deme using the low level bit
  of fwdpp
  \note In practice, one would use KTfwd::singlepop instead of this.  This
  object is for unit/integration testing only!!
  \ingroup unit
*/
{
    mcont_t mutations, fixations;
    gcont_t gametes;
    dipvector_t diploids;
    lookup_table_t mut_lookup;
    mcounts_t mcounts, fixation_times;
    KTfwd::gamete::mutation_container neutral, selected;
    gsl_rng *r;
    standard_empty_single_deme_fixture()
        : mutations(mcont_t()), fixations(mcont_t()), gametes(gcont_t()),
          diploids(dipvector_t()), mut_lookup(lookup_table_t()),
          mcounts(mcounts_t()), fixation_times(mcounts_t()),
          neutral(KTfwd::gamete::mutation_container()),
          selected(KTfwd::gamete::mutation_container()),
          r(gsl_rng_alloc(gsl_rng_mt19937))
    {
    }
    ~standard_empty_single_deme_fixture() { gsl_rng_free(r); }
};

struct standard_empty_metapop_fixture
/*!
  Basic stuff needed for a simulation of a metapopulation using the low level
  bit of fwdpp
  \note In practice, one would use KTfwd::metapop instead of this.  This object
  is for unit/integration testing only!!
  \ingroup unit
*/
{
    mcont_t mutations, fixations;
    gcont_t gametes;
    std::vector<dipvector_t> diploids;
    lookup_table_t mut_lookup;
    mcounts_t mcounts, fixation_times;
    KTfwd::gamete::mutation_container neutral, selected;
    gsl_rng *r;
    standard_empty_metapop_fixture()
        : mutations(mcont_t()), fixations(mcont_t()), gametes(gcont_t()),
          diploids(std::vector<dipvector_t>()), mut_lookup(lookup_table_t()),
          mcounts(mcounts_t()), fixation_times(mcounts_t()),
          neutral(KTfwd::gamete::mutation_container()),
          selected(KTfwd::gamete::mutation_container()),
          r(gsl_rng_alloc(gsl_rng_mt19937))
    {
    }
    ~standard_empty_metapop_fixture() { gsl_rng_free(r); }
};

struct standard_empty_multiloc_fixture
/*!
  Basic stuff needed for a simulation of a multilocus simulation using the low
  level bit of fwdpp
  \note In practice, one would use KTfwd::multilocus instead of this.  This
  object is for unit/integration testing only!!
  \ingroup unit
*/
{
    mcont_t mutations, fixations;
    gcont_t gametes;
    std::vector<dipvector_t> diploids;
    lookup_table_t mut_lookup;
    mcounts_t mcounts, fixation_times;
    KTfwd::gamete::mutation_container neutral, selected;
    gsl_rng *r;
    standard_empty_multiloc_fixture()
        : mutations(mcont_t()), fixations(mcont_t()), gametes(gcont_t()),
          diploids(std::vector<dipvector_t>()), mut_lookup(lookup_table_t()),
          mcounts(mcounts_t()), fixation_times(mcounts_t()),
          neutral(KTfwd::gamete::mutation_container()),
          selected(KTfwd::gamete::mutation_container()),
          r(gsl_rng_alloc(gsl_rng_mt19937))
    {
    }
    ~standard_empty_multiloc_fixture() { gsl_rng_free(r); }
};
#endif
