#include <vector>
#include <cstdint>
#include <stdexcept>
#include <algorithm>
#include <unordered_set>
#include <gsl/gsl_randist.h>
#include "tree_sequence_examples_types.hpp"
#include <fwdpp/ts/table_collection.hpp>
#include <fwdpp/ts/table_simplifier.hpp>
#include <fwdpp/ts/exceptions.hpp>
#include <fwdpp/ts/tree_visitor.hpp>
#include <fwdpp/ts/marginal_tree_functions.hpp>
#include <fwdpp/genetic_map/poisson_interval.hpp>

std::vector<fwdpp::ts::TS_NODE_INT>
sample_nodes_from_metadata(const std::vector<diploid_metadata>& metadata)
{
    std::vector<fwdpp::ts::TS_NODE_INT> samples;
    samples.reserve(2 * metadata.size());
    for (const auto& md : metadata)
        {
            samples.push_back(md.n1);
            samples.push_back(md.n2);
        }
    return samples;
}

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

double
get_new_mutation_position(const GSLrng& rng, double genome_length,
                          std::unordered_set<double>& mutation_positions)
{
    auto pos = gsl_ran_flat(rng.get(), 0, genome_length);
    while (mutation_positions.find(pos) != end(mutation_positions))
        {
            pos = gsl_ran_flat(rng.get(), 0, genome_length);
        }
    mutation_positions.insert(pos);
    return pos;
}

void
mutate_offspring_nodes(const GSLrng& rng, double mutrate,
                       std::size_t num_nodes_before_births,
                       std::size_t num_nodes, double genome_length,
                       fwdpp::ts::site_vector& site_table,
                       fwdpp::ts::mutation_key_vector& mutation_table,
                       std::unordered_set<double>& mutation_positions)
{
    for (std::size_t offspring_node = num_nodes_before_births;
         offspring_node < num_nodes; ++offspring_node)
        {
            unsigned n = gsl_ran_poisson(rng.get(), mutrate);
            for (std::size_t i = 0; i < n; ++i)
                {
                    auto pos = get_new_mutation_position(rng, genome_length,
                                                         mutation_positions);
                    site_table.emplace_back(fwdpp::ts::site{ pos, 0 });
                    // NOTE: fwdpp's mutation table records contain indexes
                    // to mutation metadata in population objects.  However,
                    // as these sims don't use fwdpp's mutation types, such
                    // metadata are moot.  Thus, we will record the mutation's
                    // location in the mutation table instead.
                    mutation_table.emplace_back(fwdpp::ts::mutation_record{
                        static_cast<fwdpp::ts::TS_NODE_INT>(
                            offspring_node),       // Node id
                        mutation_table.size() + 1, // See NOTE above
                        site_table.size() - 1,     // Site id
                        1,                         // derived state
                        true                       // is mutation neutral?
                    });
                }
        }
}

void
rebuld_mutation_position_lookup(const fwdpp::ts::site_vector& site_table,
                                std::unordered_set<double>& mutation_positions)
{
    mutation_positions.clear();
    for (auto&& s : site_table)
        {
            mutation_positions.insert(s.position);
        }
}

void
simplify_tables_remap_metadata(double psurvival,
                               fwdpp::ts::table_simplifier& simplifier,
                               fwdpp::ts::table_collection& tables,
                               std::vector<diploid_metadata>& metadata)
{
    // Fetch node data from metadata
    auto samples = sample_nodes_from_metadata(metadata);
    // Sort the tables
    if (psurvival > 0.0) // overlapping generations
        {
            // NOTE: this is a workaround for GitHub issue #238
            tables.edge_offset = 0;
        }
    tables.sort_tables_for_simplification();
    // Simplify
    auto x = simplifier.simplify(tables, samples);

    // Remap nodes for alive individuals
    for (auto& md : metadata)
        {
            // Validate output
            if (x.first[md.n1] == fwdpp::ts::TS_NULL_NODE
                || x.first[md.n2] == fwdpp::ts::TS_NULL_NODE)
                {
                    throw fwdpp::ts::tables_error("output node maps to null");
                }
            md.n1 = x.first[md.n1];
            md.n2 = x.first[md.n2];
        }
}

void
rebuild_edge_index(const std::vector<fwdpp::ts::indexed_edge>& current,
                   const std::vector<fwdpp::ts::indexed_edge>& new_entries,
                   std::vector<fwdpp::ts::indexed_edge>& output)
{
    output.clear();
    auto itr = begin(current);
    auto first = begin(current);
    const auto last = end(current);

    for (auto& i : new_entries)
        {
            itr = std::lower_bound(first, last, i);
            output.insert(end(output), first, itr);
            output.emplace_back(std::move(i));
            first = itr;
        }
    output.insert(end(output), first, last);
    if (output.size() != current.size() + new_entries.size())
        {
            throw std::runtime_error("rebuild_edge_index error");
        }
}

void
dynamic_update_edge_table_indexes(std::size_t num_edges,
                                  std::vector<fwdpp::ts::indexed_edge>& newI,
                                  std::vector<fwdpp::ts::indexed_edge>& newO,
                                  fwdpp::ts::table_collection& tables)
{
    if (num_edges == 0)
        {
            tables.build_indexes();
            return;
        }

    // Generate and sort the new data entries
    std::vector<fwdpp::ts::indexed_edge> I, O;
    for (auto i = begin(tables.edge_table) + num_edges;
         i < end(tables.edge_table); ++i)
        {
            I.emplace_back(i->left, -tables.node_table[i->parent].time,
                           i->parent, i->child);
            O.emplace_back(i->right, tables.node_table[i->parent].time,
                           i->parent, i->child);
        }
    std::sort(begin(I), end(I));
    std::sort(begin(O), end(O));

    rebuild_edge_index(tables.input_left, I, newI);
    rebuild_edge_index(tables.output_right, O, newO);

    tables.input_left.swap(newI);
    tables.output_right.swap(newO);
#ifndef NDEBUG
    if (!std::is_sorted(begin(tables.input_left), end(tables.input_left)))
        {
            throw std::runtime_error("input_left is not sorted");
        }
    if (!std::is_sorted(begin(tables.output_right), end(tables.output_right)))
        {
            throw std::runtime_error("output_right is not sorted");
        }
#endif
}

void
traverse_edges_check_samples(
    const fwdpp::ts::table_collection& tables,
    const std::vector<diploid_metadata>& metadata,
    std::vector<fwdpp::ts::TS_NODE_INT>& samples_buffer)
{
    auto samples = sample_nodes_from_metadata(metadata);
    fwdpp::ts::tree_visitor tv(tables, samples,
                               fwdpp::ts::update_samples_list(true));
    std::vector<int> found(tables.node_table.size(), 0);
    std::vector<int> is_sample(tables.node_table.size(), 0);
    for (auto s : samples)
        {
            is_sample[s] = 1;
        }
    const auto f = [&samples_buffer](fwdpp::ts::TS_NODE_INT u) {
        samples_buffer.push_back(u);
    };
    while (tv())
        {
            const auto& m = tv.tree();
            auto roots = fwdpp::ts::get_roots(m);
            unsigned nfound = 0;
            for (auto r : roots)
                {
                    samples_buffer.clear();
                    fwdpp::ts::process_samples(
                        m, fwdpp::ts::convert_sample_index_to_nodes(true), r,
                        f);
                    for (auto i : samples_buffer)
                        {
                            if (found[i] != 0)
                                {
                                    throw fwdpp::ts::samples_error(
                                        "sample descends from > 1 root in a "
                                        "single tree");
                                }
                            found[i]++;
                            nfound++;
                        }
                }
            if (nfound != samples.size())
                {
                    throw fwdpp::ts::samples_error(
                        "fatal error: found wrong number of samples");
                }
            // reset the data.
            for (auto i : samples)
                {
                    found[i] = 0;
                }
        }
}

void
traverse_trees_count_mutations(
    const fwdpp::ts::table_collection& tables,
    const std::vector<fwdpp::ts::TS_NODE_INT>& samples,
    const bool track_samples, std::vector<unsigned>& mutation_counts,
    std::vector<fwdpp::ts::TS_NODE_INT>& samples_buffer)
{
    std::fill(begin(mutation_counts), end(mutation_counts), 0);
    mutation_counts.resize(tables.mutation_table.size(), 0);

    fwdpp::ts::tree_visitor tv(tables, samples,
                               fwdpp::ts::update_samples_list(track_samples));
    const auto first_mutation = begin(tables.mutation_table);
    const auto last_mutation = end(tables.mutation_table);
    auto current_mutation = begin(tables.mutation_table);

    while (current_mutation < last_mutation && tv())
        {
            auto& m = tv.tree();
            // Update mutation iterator to mutations on current tree
            while (current_mutation < last_mutation
                   && tables.site_table[current_mutation->site].position
                          < m.left)
                {
                    ++current_mutation;
                }
            // Process all mutations on the current tree
            while (current_mutation < last_mutation
                   && tables.site_table[current_mutation->site].position
                          < m.right)
                {
                    mutation_counts[current_mutation - first_mutation]
                        = m.leaf_counts[current_mutation->node];
                    if (track_samples)
                        {
                            // NOTE: with dynamic indexing, there are
                            // going to be extinct mutations, so we
                            // skip those that have no descendants in the sample
                            if (m.leaf_counts[current_mutation->node])
                                {
                                    samples_buffer.clear();
                                    fwdpp::ts::process_samples(
                                        m,
                                        fwdpp::ts::
                                            convert_sample_index_to_nodes(
                                                true),
                                        current_mutation->node,
                                        [&samples_buffer](
                                            fwdpp::ts::TS_NODE_INT u) {
                                            samples_buffer.push_back(u);
                                        });
                                }
                        }
                    ++current_mutation;
                }
        }
}

void
remove_fixed_variants(unsigned twoN, fwdpp::ts::site_vector& site_table,
                      fwdpp::ts::mutation_key_vector& mutation_table,
                      std::unordered_set<double>& mutation_positions,
                      std::vector<unsigned>& mutation_counts)
{
    bool any_fixed = false;
    for (std::size_t i = 0; i < mutation_counts.size(); ++i)
        {
            if (mutation_counts[i] == twoN)
                {
                    any_fixed = true;
                    mutation_table[i].key
                        = std::numeric_limits<std::size_t>::max();
                    auto itr = mutation_positions.find(
                        site_table[mutation_table[i].site].position);
                    if (itr == end(mutation_positions))
                        {
                            throw std::runtime_error(
                                "mutation data structures are in an "
                                "inconsistent state");
                        }
                    site_table[mutation_table[i].site].position
                        = -1.; //conver to impossible value for a tree-seq sim
                }
        }
    if (any_fixed)
        {
            site_table.erase(std::remove_if(begin(site_table), end(site_table),
                                            [](fwdpp::ts::site& s) {
                                                return s.position == -1.;
                                            }),
                             end(site_table));
            mutation_table.erase(
                std::remove_if(
                    begin(mutation_table), end(mutation_table),
                    [](fwdpp::ts::mutation_record& mr) {
                        return mr.key
                               == std::numeric_limits<std::size_t>::max();
                    }),
                end(mutation_table));
            mutation_counts.erase(std::remove(begin(mutation_counts),
                                              end(mutation_counts), twoN),
                                  end(mutation_counts));
        }
}

void
wfevolvets_no_mutation(const GSLrng& rng, unsigned ngenerations,
                       unsigned simplify, double psurvival,
                       const fwdpp::poisson_interval& recombination,
                       std::vector<diploid_metadata>& metadata,
                       fwdpp::ts::table_collection& tables)
{
    if (psurvival < 0. || psurvival >= 1.0 || !std::isfinite(psurvival))
        {
            throw std::invalid_argument("invalid survival probability");
        }
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
                    simplify_tables_remap_metadata(psurvival, simplifier,
                                                   tables, metadata);
                    simplified = true;
                }
            // The current metadata become the parental metadata.
            // For efficiency, it is best if metadata are fast to copy.
            parental_metadata = metadata;
        }
    if (!simplified)
        {
            simplify_tables_remap_metadata(psurvival, simplifier, tables,
                                           metadata);
        }
}

std::vector<unsigned>
wfevolvets_dynamic_indexing(const GSLrng& rng, unsigned ngenerations,
                            unsigned count_mutations_interval,
                            bool remove_fixations, bool track_samples,
                            unsigned check_interval, unsigned simplify,
                            double psurvival, double mutrate,
                            const fwdpp::poisson_interval& recombination,
                            std::vector<diploid_metadata>& metadata,
                            fwdpp::ts::table_collection& tables)
{
    if (psurvival < 0. || psurvival >= 1.0 || !std::isfinite(psurvival))
        {
            throw std::invalid_argument("invalid survival probability");
        }
    if (mutrate < 0. || !std::isfinite(mutrate))
        {
            throw std::invalid_argument("invalid mutation rate");
        }
    fwdpp::ts::table_simplifier simplifier(tables.genome_length());

    std::vector<std::size_t> dead_individual_indexes;
    std::vector<diploid_metadata> parental_metadata(metadata);
    std::vector<double> breakpoints; // reusable buffer for rec breakpoints
    std::vector<fwdpp::ts::indexed_edge> Ibuffer, Obuffer;
    std::vector<fwdpp::ts::TS_NODE_INT> samples_buffer;
    std::unordered_set<double> mutation_positions;
    std::vector<unsigned> mutation_counts;
    bool simplified = false;
    for (unsigned gen = 0; gen < ngenerations; ++gen)
        {
            deaths(rng, metadata, psurvival, dead_individual_indexes);
            auto num_edges = tables.edge_table.size();
            // NOTE that gen + 1 will be the birth time.
            auto num_nodes_before_births = tables.node_table.size();
            births(rng, gen + 1, dead_individual_indexes, parental_metadata,
                   recombination, metadata, tables, breakpoints);
            mutate_offspring_nodes(rng, mutrate, num_nodes_before_births,
                                   tables.node_table.size(),
                                   tables.genome_length(), tables.site_table,
                                   tables.mutation_table, mutation_positions);
            dynamic_update_edge_table_indexes(num_edges, Ibuffer, Obuffer,
                                              tables);
            if ((gen + 1) % simplify == 0.0)
                {
                    simplify_tables_remap_metadata(psurvival, simplifier,
                                                   tables, metadata);
                    rebuld_mutation_position_lookup(tables.site_table,
                                                    mutation_positions);
                    tables.build_indexes();
                    simplified = true;
                }
            if (check_interval && (gen + 1) % check_interval == 0.0)
                {
                    traverse_edges_check_samples(tables, metadata,
                                                 samples_buffer);
                }
            if (count_mutations_interval
                && (gen + 1) % count_mutations_interval == 0.0)
                {
                    auto samples = sample_nodes_from_metadata(metadata);
                    traverse_trees_count_mutations(
                        tables, samples, track_samples, mutation_counts,
                        samples_buffer);
                    if (remove_fixations)
                        {
                            remove_fixed_variants(
                                2 * metadata.size(), tables.site_table,
                                tables.mutation_table, mutation_positions,
                                mutation_counts);
                        }
                }
            // The current metadata become the parental metadata.
            // For efficiency, it is best if metadata are fast to copy.
            parental_metadata = metadata;
        }
    if (!simplified)
        {
            simplify_tables_remap_metadata(psurvival, simplifier, tables,
                                           metadata);
            rebuld_mutation_position_lookup(tables.site_table,
                                            mutation_positions);
            tables.build_indexes();
        }
    if (count_mutations_interval)
        {

            auto samples = sample_nodes_from_metadata(metadata);
            traverse_trees_count_mutations(tables, samples, track_samples,
                                           mutation_counts, samples_buffer);
            if (remove_fixations)
                {
                    remove_fixed_variants(2 * metadata.size(),
                                          tables.site_table,
                                          tables.mutation_table,
                                          mutation_positions, mutation_counts);
                }
            // There can be no mutations in the table with count == 0
            for (std::size_t i = 0; i < tables.mutation_table.size(); ++i)
                {
                    if (mutation_counts[i] == 0)
                        {
                            throw std::runtime_error(
                                "mutation count cannot be zero");
                        }
                    if (mutation_counts[i] > 2 * metadata.size())
                        {
                            throw std::runtime_error(
                                "mutation count is greater than number of "
                                "genomes");
                        }
                }
        }
    return mutation_counts;
}
