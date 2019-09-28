#include <vector>
#include <cstdint>
#include <stdexcept>
#include <algorithm>
#include <gsl/gsl_randist.h>
#include "tree_sequence_examples_types.hpp"
#include <fwdpp/ts/table_collection.hpp>
#include <fwdpp/ts/table_simplifier.hpp>
#include <fwdpp/ts/exceptions.hpp>
#include <fwdpp/genetic_map/poisson_interval.hpp>

void
deaths(const GSLrng& rng, const std::vector<diploid_metadata>& metadata,
       double psurvival, std::vector<std::size_t>& dead_individual_indexes)
{
    dead_individual_indexes.clear();
    for (std::size_t i = 0; i < metadata.size(); ++i)
        {
            if (psurvival == 0.0 || gsl_rng_uniform(rng.get()) > psurvival)
                {
                    dead_individual_indexes.push_back(i);
                }
        }
}

std::tuple<fwdpp::ts::TS_NODE_INT, fwdpp::ts::TS_NODE_INT>
pick_parent(const GSLrng& rng,
            const std::vector<diploid_metadata>& parental_metadata)
{
    std::size_t p = gsl_ran_flat(rng.get(), 0, parental_metadata.size());
    auto n1 = parental_metadata[p].n1;
    auto n2 = parental_metadata[p].n2;

    if (gsl_rng_uniform(rng.get()) < 0.5)
        {
            std::swap(n1, n2);
        }
    return std::make_tuple(n1, n2);
}

void
generate_recombination_breakpoints(
    const GSLrng& rng, const fwdpp::poisson_interval& recombination,
    std::vector<double>& breakpoints)
{
    breakpoints.clear();
    recombination(rng.get(), breakpoints);
    if (!breakpoints.empty())
        {
            std::sort(begin(breakpoints), end(breakpoints));
            breakpoints.push_back(std::numeric_limits<double>::max());
        }
}

void
births(const GSLrng& rng, const unsigned birth_time,
       const std::vector<std::size_t>& dead_individual_indexes,
       const std::vector<diploid_metadata>& parental_metadata,
       const fwdpp::poisson_interval& recombination,
       std::vector<diploid_metadata>& metadata,
       fwdpp::ts::table_collection& tables, std::vector<double>& breakpoints)
{
    for (std::size_t b = 0; b < dead_individual_indexes.size(); ++b)
        {
            auto parental_nodes = pick_parent(rng, parental_metadata);
            generate_recombination_breakpoints(rng, recombination,
                                               breakpoints);
            auto offspring_first_node = tables.register_diploid_offspring(
                breakpoints, parental_nodes, 0, birth_time);
            parental_nodes = pick_parent(rng, parental_metadata);
            generate_recombination_breakpoints(rng, recombination,
                                               breakpoints);
            auto offspring_second_node = tables.register_diploid_offspring(
                breakpoints, parental_nodes, 0, birth_time);

            // Take a reference to the metadata for the individual
            // we are replacing.  This reference streamlines the
            // code for updating metadata
            auto& md = metadata[dead_individual_indexes[b]];
            md.time = birth_time;
            md.n1 = offspring_first_node;
            md.n2 = offspring_second_node;
        }
}

void
simplify_tables_remap_metadata(fwdpp::ts::table_simplifier& simplifier,
                               fwdpp::ts::table_collection& tables,
                               std::vector<diploid_metadata>& metadata)
{
    // Fetch node data from metadata
    std::vector<fwdpp::ts::TS_NODE_INT> samples;
    samples.reserve(2 * metadata.size());
    for (const auto& md : metadata)
        {
            samples.push_back(md.n1);
            samples.push_back(md.n2);
        }
    // Sort the tables
    tables.sort_tables_for_simplification();
    // Simplify
    auto x = simplifier.simplify(tables, samples);

    // Remap nodes for alive individuals
    for (auto& md : metadata)
        {
            // Validate output
            if (x.first[md.n1] == fwdpp::ts::TS_NULL_NODE
                || x.first[md.n1] == fwdpp::ts::TS_NULL_NODE)
                {
                    throw fwdpp::ts::tables_error("output node maps to null");
                }
            md.n1 = x.first[md.n1];
            md.n2 = x.first[md.n2];
        }
}

void
wfevolvets_no_mutation(const GSLrng& rng, unsigned ngenerations,
                       unsigned simplify, double psurvival,
                       const fwdpp::poisson_interval& recombination,
                       std::vector<diploid_metadata>& metadata,
                       fwdpp::ts::table_collection& tables)
{
    // TODO: validate input params
    fwdpp::ts::table_simplifier simplifier(tables.genome_length());

    std::vector<std::size_t> dead_individual_indexes;
    std::vector<diploid_metadata> parental_metadata(metadata);
    std::vector<double> breakpoints; // reusable buffer for rec breakpoints
    bool simplified = false;
    for (unsigned gen = 0; gen < ngenerations; ++gen)
        {
            simplified = false;
            deaths(rng, metadata, psurvival, dead_individual_indexes);
            // NOTE that gen + 1 will be the birth time.
            births(rng, gen + 1, dead_individual_indexes, parental_metadata,
                   recombination, metadata, tables, breakpoints);
            if ((gen + 1) % simplify == 0.0)
                {
                    simplify_tables_remap_metadata(simplifier, tables,
                                                   metadata);
                    simplified = true;
                }
            // The current metadata become the parental metadata.
            // For efficiency, it is best if metadata are fast to copy.
            parental_metadata = metadata;
        }
    if (!simplified)
        {
            simplify_tables_remap_metadata(simplifier, tables, metadata);
        }
}
