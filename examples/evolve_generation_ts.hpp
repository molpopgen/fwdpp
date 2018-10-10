// Wow, that's a lot of stuff needed:
template <typename breakpoint_function, typename mutation_model,
          typename mrecbin, typename grecbin>
std::int32_t
generate_offspring(const GSLrng_t& rng, const breakpoint_function& recmodel,
                   const mutation_model& mmodel, const double mu,
                   const std::size_t parent, const fwdpp::uint_t parent_g1,
                   const fwdpp::uint_t parent_g2,
                   const std::tuple<std::int32_t, std::int32_t>& parent_nodes,
                   const std::int32_t generation,
                   const std::int32_t next_index, slocuspop_t& pop,
                   std::size_t& offspring_gamete, table_collection& tables,
                   mrecbin& mutation_recycling_bin,
                   grecbin& gamete_recycling_bin)
{
    auto breakpoints = recmodel();
    auto new_mutations = fwdpp::generate_new_mutations(
        mutation_recycling_bin, rng.get(), mu, pop.diploids[parent],
        pop.gametes, pop.mutations, parent_g1, mmodel);
#ifndef NDEBUG
    for (auto& m : new_mutations)
        {
            auto itr = pop.mut_lookup.equal_range(pop.mutations[m].pos);
            assert(std::distance(itr.first, itr.second) == 1);
        }
#endif
    // We will only add selected mutations into offspring gametes.
    auto end_of_neutral
        = std::stable_partition(new_mutations.begin(), new_mutations.end(),
                                [&pop](const fwdpp::uint_t key) {
                                    return pop.mutations[key].neutral == true;
                                });
    //if (!std::is_sorted(end_of_neutral, new_mutations.end(),
    //                    [&pop](const std::size_t a, const std::size_t b) {
    //                        return pop.mutations[a].pos < pop.mutations[b].pos;
    //                    }))
    //    {
    //        throw std::runtime_error("bad");
    //    }
    offspring_gamete = fwdpp::mutate_recombine(
        decltype(new_mutations)(end_of_neutral, new_mutations.end()),
        breakpoints, parent_g1, parent_g2, pop.gametes, pop.mutations,
        gamete_recycling_bin, pop.neutral, pop.selected);
    //if (!new_mutations.empty() || !breakpoints.empty())
    //    {
    //        assert(offspring_gamete != parent_g1);
    //    }
    tables.add_offspring_data(next_index, breakpoints, new_mutations,
                              parent_nodes, generation);
    return next_index + 1;
}

template <typename breakpoint_function, typename mutation_model>
void
evolve_generation(const GSLrng_t& rng, slocuspop_t& pop,
                  const fwdpp::uint_t N_next, const double mu,
                  const mutation_model& mmodel,
                  const breakpoint_function& recmodel,
                  const fwdpp::uint_t generation, table_collection& tables,
                  table_simplifier& simplifier,
                  std::int32_t first_parental_index, std::int32_t next_index)
{

    auto gamete_recycling_bin
        = fwdpp::fwdpp_internal::make_gamete_queue(pop.gametes);
    auto mutation_recycling_bin
        = fwdpp::fwdpp_internal::make_mut_queue(pop.mcounts);

    auto lookup = w(pop, fwdpp::multiplicative_diploid());
    decltype(pop.diploids) offspring(N_next);

    // Generate the offspring
    auto next_index_local = next_index;
    for (auto& dip : offspring)
        {
            auto p1 = gsl_ran_discrete(rng.get(), lookup.get());
            auto p2 = gsl_ran_discrete(rng.get(), lookup.get());
            auto p1g1 = pop.diploids[p1].first;
            auto p1g2 = pop.diploids[p1].second;
            auto p2g1 = pop.diploids[p2].first;
            auto p2g2 = pop.diploids[p2].second;

            // Mendel
            int swap1 = (gsl_rng_uniform(rng.get()) < 0.5) ? 1 : 0;
            int swap2 = (gsl_rng_uniform(rng.get()) < 0.5) ? 1 : 0;
            if (swap1)
                std::swap(p1g1, p1g2);
            if (swap2)
                std::swap(p2g1, p2g2);

            auto p1id = get_parent_ids(first_parental_index, p1, swap1);
            auto p2id = get_parent_ids(first_parental_index, p2, swap2);

            assert(std::get<0>(p1id) < 2 * static_cast<std::int32_t>(N_next));
            assert(std::get<1>(p1id) < 2 * static_cast<std::int32_t>(N_next));
            assert(std::get<0>(p2id) < 2 * static_cast<std::int32_t>(N_next));
            assert(std::get<1>(p2id) < 2 * static_cast<std::int32_t>(N_next));

            next_index_local = generate_offspring(
                rng, recmodel, mmodel, mu, p1, p1g1, p1g2, p1id, generation,
                next_index_local, pop, dip.first, tables,
                mutation_recycling_bin, gamete_recycling_bin);
            next_index_local = generate_offspring(
                rng, recmodel, mmodel, mu, p2, p2g1, p2g2, p2id, generation,
                next_index_local, pop, dip.second, tables,
                mutation_recycling_bin, gamete_recycling_bin);
            pop.gametes[dip.first].n++;
            pop.gametes[dip.second].n++;
        }
    assert(next_index_local
           == next_index + 2 * static_cast<std::int32_t>(N_next));
    // This is constant-time
    pop.diploids.swap(offspring);
    tables.sort_tables(pop.mutations);
    std::vector<std::int32_t> samples(2 * pop.diploids.size());
    std::iota(samples.begin(), samples.end(),
              tables.num_nodes() - 2 * pop.diploids.size());
    auto idmap = simplifier.simplify(tables, samples, pop.mutations);
    tables.build_indexes();
    for (auto& s : samples)
        {
            s = idmap[s];
        }
    tables.count_mutations(pop.mutations, samples, pop.mcounts);
#ifndef NDEBUG
    std::vector<std::size_t> keys;
    for (auto& mr : tables.mutation_table)
        {
            keys.push_back(mr.key);
        }
    std::sort(keys.begin(), keys.end());
    auto u = std::unique(keys.begin(), keys.end());
    ;
    if (u != keys.end())
        {
            std::cout << "redundant keys " << generation << '\n';
        }

    for (auto& mr : tables.mutation_table)
        {
            assert(mr.node != fwdpp::ts::TS_NULL_NODE);
            assert(tables.node_table[mr.node].generation
                   >= pop.mutations[mr.key].g);
        }
    decltype(pop.mcounts) mc;
    fwdpp::fwdpp_internal::process_gametes(pop.gametes, pop.mutations, mc);
    //assert(pop.mcounts == mc);
    for (std::size_t i = 0; i < pop.mcounts.size(); ++i)
        {
            if (!pop.mutations[i].neutral)
                {
                    assert(pop.mcounts[i] == mc[i]);
                }
        }
#endif
    tables.mutation_table.erase(
        std::remove_if(
            tables.mutation_table.begin(), tables.mutation_table.end(),
            [&pop](const fwdpp::ts::mutation_record& mr) {
                return pop.mcounts[mr.key] == 2 * pop.diploids.size();
            }),
        tables.mutation_table.end());
    fwdpp::fwdpp_internal::gamete_cleaner(
        pop.gametes, pop.mutations, pop.mcounts, 2 * N_next, std::true_type());
    fwdpp::update_mutations(pop.mutations, pop.fixations, pop.fixation_times,
                            pop.mut_lookup, pop.mcounts, generation,
                            2 * pop.diploids.size());
    if (generation && generation % 100 == 0.0)
        {
            auto new_mut_indexes = fwdpp::compact_mutations(pop);
            // Mutation compacting re-orders the data, so we need to
            // re-index our mutation table
            for (auto& mr : tables.mutation_table)
                {
                    mr.key = new_mut_indexes[mr.key];
                }
        }
}

