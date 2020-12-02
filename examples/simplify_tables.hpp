#ifndef FWDPP_EXAMPLES_SIMPLIFY_TABLES_HPP
#define FWDPP_EXAMPLES_SIMPLIFY_TABLES_HPP

#include <cstdint>
#include <vector>
#include <fwdpp/ts/std_table_collection.hpp>
#include <fwdpp/ts/table_collection_functions.hpp>
#include <fwdpp/ts/table_simplifier.hpp>
#include <fwdpp/ts/count_mutations.hpp>
#include <fwdpp/ts/remove_fixations_from_gametes.hpp>
#include <fwdpp/ts/recycling.hpp>
#include "confirm_mutation_counts.hpp"
#include "fwdpp/ts/definitions.hpp"
#include "fwdpp/ts/exceptions.hpp"

template <typename Simplifier, typename poptype>
std::pair<std::vector<fwdpp::ts::table_index_t>, std::vector<std::size_t>>
simplify_tables(poptype &pop, const fwdpp::uint_t generation,
                std::vector<fwdpp::uint_t> &mcounts_from_preserved_nodes,
                fwdpp::ts::std_table_collection &tables, Simplifier &simplifier,
                const fwdpp::ts::table_index_t first_sample_node,
                const std::size_t num_samples,
                std::vector<fwdpp::ts::table_index_t> &preserved_nodes,
                const bool preserve_fixations = false)
{
    fwdpp::ts::sort_tables_for_simplification(tables.edge_offset, tables);
    std::vector<std::int32_t> samples(num_samples);
    std::iota(samples.begin(), samples.end(), first_sample_node);
    samples.insert(end(samples), begin(preserved_nodes), end(preserved_nodes));
    auto rv = simplifier.simplify(tables, samples);
    tables.edge_offset = tables.num_edges();
    tables.build_indexes();
    for (auto &s : samples)
        {
            s = rv.first[s];
            assert(s != fwdpp::ts::TS_NULL_NODE);
        }
    for (auto &s : preserved_nodes)
        {
            s = rv.first[s];
            assert(s != fwdpp::ts::TS_NULL_NODE);
        }
    if (!preserve_fixations)
        {
            samples.resize(num_samples);
            fwdpp::ts::count_mutations(tables, pop.mutations, samples, preserved_nodes,
                                       pop.mcounts, mcounts_from_preserved_nodes);
            auto itr = std::remove_if(
                tables.mutations.begin(), tables.mutations.end(),
                [&pop,
                 &mcounts_from_preserved_nodes](const fwdpp::ts::mutation_record &mr) {
                    return pop.mcounts[mr.key] == 2 * pop.diploids.size()
                           && mcounts_from_preserved_nodes[mr.key] == 0;
                });
            auto d = std::distance(itr, end(tables.mutations));
            tables.mutations.erase(itr, tables.mutations.end());
            if (d)
                {
                    fwdpp::ts::rebuild_site_table(tables);
                }
            confirm_mutation_counts(pop, tables);
            fwdpp::ts::remove_fixations_from_haploid_genomes(
                pop.haploid_genomes, pop.mutations, pop.mcounts,
                mcounts_from_preserved_nodes, 2 * pop.diploids.size(), false);

            fwdpp::ts::flag_mutations_for_recycling(
                pop, mcounts_from_preserved_nodes, 2 * pop.diploids.size(), generation,
                std::false_type(), std::false_type());
            confirm_mutation_counts(pop, tables);
        }
    else // Need to remove non-preserved variants from the hash table
        {
            // NOTE: simplification returns the preserved node indexes
            // and does not know anything about the total number of mutations,
            // so we use a linear time/linear memory method to removed
            // extinct variants from the simulation.
            std::vector<int> preserved(pop.mutations.size(), 0);
            for (auto x : rv.second)
                {
                    preserved[x] = 1;
                }
            for (std::size_t p = 0; p < preserved.size(); ++p)
                {
                    if (!preserved[p])
                        {
                            fwdpp::ts::detail::process_mutation_index(pop.mutations,
                                                                      pop.mut_lookup, p);
                        }
                }
        }

    return rv;
}
#endif
