#ifndef FWDPP_EXAMPLES_WFEVOLVETS_HPP
#define FWDPP_EXAMPLES_WFEVOLVETS_HPP

#include "tree_sequence_examples_types.hpp"
#include <fwdpp/ts/table_collection.hpp>
#include <fwdpp/genetic_map/poisson_interval.hpp>

void wfevolvets_no_mutation(const GSLrng& rng, unsigned ngenerations,
                            unsigned simplify, double psurvival,
                            const fwdpp::poisson_interval& recombination,
                            std::vector<diploid_metadata>& metadata,
                            fwdpp::ts::table_collection& tables);

std::vector<unsigned>
wfevolvets_dynamic_indexing(const GSLrng& rng, unsigned ngenerations,
                            unsigned count_mutations_interval, bool remove_fixations,
                            bool track_samples,
                            unsigned check_interval, unsigned simplify,
                            double psurvival, double mutrate,
                            const fwdpp::poisson_interval& recombination,
                            std::vector<diploid_metadata>& metadata,
                            fwdpp::ts::table_collection& tables);

#endif

