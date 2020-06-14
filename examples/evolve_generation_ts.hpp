#ifndef FWDPP_EXAMPLES_EVOLVE_GENERATION_TS
#define FWDPP_EXAMPLES_EVOLVE_GENERATION_TS

#include <cstdint>
#include <algorithm>
#include <vector>
#include <tuple>
#include <gsl/gsl_randist.h>

#include <fwdpp/util/wrapped_range.hpp>
#include <fwdpp/util.hpp>
#include <fwdpp/internal/haploid_genome_cleaner.hpp>
#include <fwdpp/mutate_recombine.hpp>
#include <fwdpp/ts/get_parent_ids.hpp>
#include <fwdpp/ts/table_collection.hpp>
#include <fwdpp/ts/table_simplifier.hpp>
#include <fwdpp/ts/generate_offspring.hpp>
#include <fwdpp/ts/recording/mutations.hpp>
#include <fwdpp/ts/recording/diploid_offspring.hpp>

template <typename poptype, typename rng_t, typename genetic_param_holder>
std::pair<fwdpp::ts::mut_rec_intermediates, fwdpp::ts::mut_rec_intermediates>
generate_offspring(const rng_t& rng,
                   const std::pair<std::size_t, std::size_t> parent_indexes,
                   poptype& pop, typename poptype::diploid_t& offspring,
                   genetic_param_holder& genetics)
{
    auto offspring_data = fwdpp::ts::generate_offspring(
        rng.get(), parent_indexes, fwdpp::ts::selected_variants_only(), pop,
        genetics, offspring);
#ifndef NDEBUG
    for (auto& m : offspring_data.first.mutation_keys)
        {
            auto itr = pop.mut_lookup.equal_range(pop.mutations[m].pos);
            assert(std::distance(itr.first, itr.second) == 1);
        }
    for (auto& m : offspring_data.second.mutation_keys)
        {
            auto itr = pop.mut_lookup.equal_range(pop.mutations[m].pos);
            assert(std::distance(itr.first, itr.second) == 1);
        }
#endif
    return offspring_data;
}

template <typename rng_t, typename poptype, typename pick_parent1_fxn,
          typename pick_parent2_fxn, typename offspring_metadata_fxn,
          typename genetic_param_holder>
void
evolve_generation(const rng_t& rng, poptype& pop,
                  genetic_param_holder& genetics, const fwdpp::uint_t N_next,
                  const pick_parent1_fxn& pick1, const pick_parent2_fxn& pick2,
                  const offspring_metadata_fxn& update_offspring,
                  const fwdpp::uint_t generation,
                  fwdpp::ts::std_table_collection& tables,
                  std::int32_t first_parental_index, std::int32_t next_index)
{
    fwdpp::debug::all_haploid_genomes_extant(pop);

    genetics.haploid_genome_recycling_bin
        = fwdpp::make_haploid_genome_queue(pop.haploid_genomes);

    fwdpp::zero_out_haploid_genomes(pop);

    decltype(pop.diploids) offspring(N_next);

    // Generate the offspring
    auto next_index_local = next_index;
    for (std::size_t next_offspring = 0; next_offspring < offspring.size();
         ++next_offspring)
        {
            auto p1 = pick1();
            auto p2 = pick2(p1);
            auto& dip = offspring[next_offspring];
            auto offspring_data = generate_offspring(
                rng, std::make_pair(p1, p2), pop, dip, genetics);
            auto p1id = fwdpp::ts::get_parent_ids(
                first_parental_index, p1, offspring_data.first.swapped);
            auto p2id = fwdpp::ts::get_parent_ids(
                first_parental_index, p2, offspring_data.second.swapped);
            next_index_local = fwdpp::ts::record_diploid_offspring(
                offspring_data.first.breakpoints, p1id, 0, generation, tables);
            fwdpp::ts::record_mutations_infinite_sites(
                next_index_local, pop.mutations,
                offspring_data.first.mutation_keys, tables);
            next_index_local = fwdpp::ts::record_diploid_offspring(
                offspring_data.second.breakpoints, p2id, 0, generation, tables);
            fwdpp::ts::record_mutations_infinite_sites(
                next_index_local, pop.mutations,
                offspring_data.second.mutation_keys, tables);

            // Give the caller a chance to generate
            // any metadata for the offspring that
            // may depend on the parents
            update_offspring(next_offspring, p1, p2);
        }
    assert(next_index_local
           == next_index + 2 * static_cast<std::int32_t>(N_next) - 1);
    // This is constant-time
    pop.diploids.swap(offspring);
}

#endif
