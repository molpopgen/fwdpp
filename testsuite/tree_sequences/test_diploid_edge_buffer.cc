#include <numeric>
#include <boost/test/unit_test.hpp>
#include <gsl/gsl_randist.h>
#include <fwdpp/GSLrng_t.hpp>
#include <fwdpp/ts/diploid_edge_buffer.hpp>
#include <fwdpp/ts/table_simplifier.hpp>
#include <iostream>

std::pair<fwdpp::ts::TS_NODE_INT, unsigned>
diploid_wf(const fwdpp::GSLrng_mt& rng, unsigned N,
           unsigned parental_generation, unsigned ngenerations,
           fwdpp::ts::TS_NODE_INT first_parental_index,
           fwdpp::ts::table_collection& tables,
           fwdpp::ts::diploid_edge_buffer& buffer, bool buffer_edges)
// Based on the tskit tutorial examples https://tskit-dev.github.io/tutorials/wfforward.html
{
    decltype(first_parental_index) parent_index_offset = first_parental_index;
    fwdpp::ts::TS_NODE_INT next_offspring_node;
    std::vector<double> breakpoints; //No rec.
    const fwdpp::ts::TS_NODE_INT offspring_deme{ 0 };
    unsigned generation = 1;
    for (; generation <= ngenerations; ++generation)
        {
            if (buffer_edges)
                {
                    buffer.begin_epoch(N);
                }
            for (unsigned offspring = 0; offspring < N; ++offspring)
                {
                    std::size_t parent = gsl_ran_flat(rng.get(), 0, N);
                    fwdpp::ts::TS_NODE_INT pn0
                        = parent_index_offset + 2 * parent;
                    fwdpp::ts::TS_NODE_INT pn1 = pn0 + 1;
                    if (gsl_rng_uniform(rng.get()) < 0.5)
                        {
                            std::swap(pn0, pn1);
                        }
                    if (buffer_edges)
                        {
                            next_offspring_node = fwdpp::ts::
                                register_offspring_and_buffer_edges(
                                    breakpoints, parent,
                                    std::make_tuple(pn0, pn1), offspring_deme,
                                    parental_generation + generation, tables,
                                    buffer);
                        }
                    else
                        {
                            next_offspring_node
                                = tables.register_diploid_offspring(
                                    breakpoints, std::make_tuple(pn0, pn1),
                                    offspring_deme,
                                    parental_generation + generation);
                        }
                    if (next_offspring_node == fwdpp::ts::TS_NULL_NODE)
                        {
                            // NOTE: this is impossible, and really just
                            // silences unused variable warnings re:
                            // next_offspring_node
                            throw std::runtime_error("NULL offspring node");
                        }

                    parent = gsl_ran_flat(rng.get(), 0, N);
                    pn0 = parent_index_offset + 2 * parent;
                    pn1 = pn0 + 1;
                    if (gsl_rng_uniform(rng.get()) < 0.5)
                        {
                            std::swap(pn0, pn1);
                        }

                    if (buffer_edges)
                        {
                            next_offspring_node = fwdpp::ts::
                                register_offspring_and_buffer_edges(
                                    breakpoints, parent,
                                    std::make_tuple(pn0, pn1), offspring_deme,
                                    parental_generation + generation, tables,
                                    buffer);
                        }
                    else
                        {
                            next_offspring_node
                                = tables.register_diploid_offspring(
                                    breakpoints, std::make_tuple(pn0, pn1),
                                    offspring_deme,
                                    parental_generation + generation);
                        }
                    if (next_offspring_node == fwdpp::ts::TS_NULL_NODE)
                        {
                            throw std::runtime_error("NULL offspring node");
                        }
                }
            parent_index_offset += 2 * N;
            if (buffer_edges)
                {
                    buffer.end_epoch();
                }
        }
    return std::make_pair(parent_index_offset,
                          parental_generation + generation - 1);
}

using parent_nodes
    = std::tuple<fwdpp::ts::TS_NODE_INT, fwdpp::ts::TS_NODE_INT>;

struct parent_record
{
    std::size_t index;
    parent_nodes nodes;
    parent_record(std::size_t i, fwdpp::ts::TS_NODE_INT a,
                  fwdpp::ts::TS_NODE_INT b)
        : index(i), nodes{ a, b }
    {
    }
};

void
simple_overlapping_generations(const fwdpp::GSLrng_mt& rng, unsigned K,
                               unsigned nbirths_per_step, unsigned nsteps,
                               std::vector<parent_record>& alive_parents,
                               fwdpp::ts::table_collection& tables,
                               fwdpp::ts::diploid_edge_buffer& buffer)
// Simplistic simulation of overlapping generations.
// While this works, it reveals deficiencies in the edge_buffer
// interface.  The notes below mark where the API could provide
// better abstractions
{
    buffer.begin_epoch(alive_parents.size());

    auto next_parental_index = alive_parents.size();
    std::vector<double> breakpoints; //no recombination
    std::vector<std::size_t> alive, will_die;
    for (unsigned step = 0; step < nsteps; ++step)
        {
            std::size_t last_parent = alive_parents.size();
            std::size_t pre_birth_epoch_size = buffer.current_epoch.size();
            for (unsigned birth = 0; birth < nbirths_per_step; ++birth)
                {
                    std::size_t parent
                        = gsl_ran_flat(rng.get(), 0, last_parent);
                    if (alive_parents[parent].index
                        == std::numeric_limits<std::size_t>::max())
                        {
                            throw std::runtime_error("dead parent!");
                        }
                    auto pn0 = std::get<0>(alive_parents[parent].nodes);
                    auto pn1 = std::get<1>(alive_parents[parent].nodes);
                    if (gsl_rng_uniform(rng.get()) < 0.5)
                        {
                            std::swap(pn0, pn1);
                        }
                    auto next_offspring_node1
                        = fwdpp::ts::register_offspring_and_buffer_edges(
                            breakpoints, alive_parents[parent].index,
                            std::make_tuple(pn0, pn1), 0, step + 1, tables,
                            buffer);
                    parent = gsl_ran_flat(rng.get(), 0, last_parent);
                    if (alive_parents[parent].index
                        == std::numeric_limits<std::size_t>::max())
                        {
                            throw std::runtime_error("dead parent!");
                        }
                    pn0 = std::get<0>(alive_parents[parent].nodes);
                    pn1 = std::get<1>(alive_parents[parent].nodes);
                    if (gsl_rng_uniform(rng.get()) < 0.5)
                        {
                            std::swap(pn0, pn1);
                        }
                    auto next_offspring_node2
                        = fwdpp::ts::register_offspring_and_buffer_edges(
                            breakpoints, alive_parents[parent].index,
                            std::make_tuple(pn0, pn1), 0, step + 1, tables,
                            buffer);
                    alive_parents.emplace_back(next_parental_index++,
                                               next_offspring_node1,
                                               next_offspring_node2);
                    // NOTE: poor abstraction
                    buffer.current_epoch.resize(buffer.current_epoch.size()
                                                + 1);
                }
            // NOTE: poor abstraction
            buffer.offsets.emplace_back(pre_birth_epoch_size,
                                        buffer.current_epoch.size());
            // Now, regulate
            alive.resize(alive_parents.size());
            will_die.resize(alive_parents.size() - K);
            std::iota(begin(alive), end(alive), 0);
            gsl_ran_choose(rng.get(), will_die.data(), will_die.size(),
                           alive.data(), alive.size(), sizeof(std::size_t));
            for (auto wd : will_die)
                {
                    alive_parents[wd].index
                        = std::numeric_limits<std::size_t>::max();
                }
            alive_parents.erase(
                std::remove_if(
                    begin(alive_parents), end(alive_parents),
                    [](const parent_record& pr) {
                        return pr.index
                               == std::numeric_limits<std::size_t>::max();
                    }),
                end(alive_parents));
        }
    // NOTE: poor abstraction
    for (auto i = buffer.offsets.rbegin(); i != buffer.offsets.rend(); ++i)
        {
            for (auto j = begin(buffer.current_epoch) + i->first;
                 j < begin(buffer.current_epoch) + i->second; ++j)
                {
                    for (auto&& k : *j)
                        {
                            tables.edge_table.insert(end(tables.edge_table),
                                                     begin(k), end(k));
                            k.clear();
                        }
                }
        }
}

struct edge_buffer_fixture
{
    const double L;
    const fwdpp::ts::TS_NODE_INT N;
    fwdpp::ts::table_collection tables, tables2;
    fwdpp::ts::diploid_edge_buffer buffer, buffer2;
    // NOTE: have to have two rng
    // when comparing output from
    // stochastic functions.
    fwdpp::GSLrng_mt rng, rng2;
    unsigned parental_generation;
    edge_buffer_fixture()
        : L{ 1.0 }, N{ 1000 }, tables(2 * N, 0.0, 0, L),
          tables2(tables), buffer{}, buffer2{}, rng{ 42 }, rng2{ 42 },
          parental_generation{ 0 }
    {
    }
};

class edge_buffer_fixture_for_overlapping : public edge_buffer_fixture
{
  private:
    std::vector<parent_record>
    fill_parents(fwdpp::ts::TS_NODE_INT num_individuals)
    {
        std::vector<parent_record> rv;
        for (decltype(num_individuals) i = 0; i < num_individuals; ++i)
            {
                rv.emplace_back(i, 2 * i, 2 * i + 1);
            }
        return rv;
    }

  public:
    std::vector<parent_record> alive_parents;
    edge_buffer_fixture_for_overlapping()
        : edge_buffer_fixture(), alive_parents(fill_parents(N))
    {
    }
};

BOOST_AUTO_TEST_SUITE(test_diploid_edge_buffer)

BOOST_FIXTURE_TEST_CASE(test_sort_order, edge_buffer_fixture)
{
    auto p = diploid_wf(rng, N, parental_generation, 100, 0, tables, buffer,
                        true);
    parental_generation = p.second;
    buffer.prepare_tables_for_simplification(tables);
    BOOST_REQUIRE_EQUAL(tables.edge_table.size(), parental_generation * 2 * N);
    BOOST_REQUIRE_EQUAL(tables.node_table.size(),
                        (parental_generation + 1) * 2 * N);
    BOOST_REQUIRE_EQUAL(tables.edges_are_sorted(), true);
}

BOOST_FIXTURE_TEST_CASE(test_two_rounds, edge_buffer_fixture)
{
    auto p = diploid_wf(rng, N, parental_generation, 100, 0, tables, buffer,
                        true);
    buffer.prepare_tables_for_simplification(tables);
    tables.edge_offset = tables.edge_table.size();
    parental_generation = p.second;
    p = diploid_wf(rng, N, parental_generation, 100, p.first, tables, buffer,
                   true);
    parental_generation = p.second;
    buffer.prepare_tables_for_simplification(tables);
    BOOST_REQUIRE_EQUAL(tables.edge_table.size(), parental_generation * 2 * N);
    BOOST_REQUIRE_EQUAL(tables.node_table.size(),
                        (parental_generation + 1) * 2 * N);
    BOOST_REQUIRE_EQUAL(tables.edges_are_sorted(), true);
}

BOOST_FIXTURE_TEST_CASE(test_simplification, edge_buffer_fixture)
{
    auto p = diploid_wf(rng, N, parental_generation, 100, 0, tables, buffer,
                        true);
    buffer.prepare_tables_for_simplification(tables);
    BOOST_CHECK_EQUAL(tables.edges_are_sorted(), true);
    // NOTE: have to use the other rng here.
    auto p2 = diploid_wf(rng2, N, parental_generation, 100, 0, tables2,
                         buffer2, false);
    tables2.sort_tables();
    BOOST_CHECK(p == p2);
    BOOST_CHECK_EQUAL(tables2.edges_are_sorted(), true);
    BOOST_REQUIRE(tables == tables2);

    std::vector<fwdpp::ts::TS_NODE_INT> samples(2 * N);
    std::iota(begin(samples), end(samples), p.first);

    fwdpp::ts::table_simplifier simplifier(L);

    auto output = simplifier.simplify(tables, samples);
    auto output_copy = simplifier.simplify(tables2, samples);
    BOOST_CHECK_EQUAL(tables2.edges_are_minimally_sorted(), true);
    BOOST_CHECK(output == output_copy);
    BOOST_CHECK(tables == tables2);
}

BOOST_FIXTURE_TEST_CASE(test_overlapping_generations,
                        edge_buffer_fixture_for_overlapping)
{
    simple_overlapping_generations(rng, N, 100, 100, alive_parents, tables,
                                   buffer);
    BOOST_CHECK(tables.edges_are_sorted());
}

BOOST_FIXTURE_TEST_CASE(test_overlapping_generations_simplify,
                        edge_buffer_fixture_for_overlapping)
{
    simple_overlapping_generations(rng, N, 100, 100, alive_parents, tables,
                                   buffer);
    BOOST_CHECK(tables.edges_are_sorted());
    std::vector<fwdpp::ts::TS_NODE_INT> samples;
    for (auto&& a : alive_parents)
        {
            samples.push_back(std::get<0>(a.nodes));
            samples.push_back(std::get<1>(a.nodes));
        }
    fwdpp::ts::table_simplifier simplifier(tables.genome_length());
    auto rv = simplifier.simplify(tables, samples);
    for (auto s : samples)
        {
            BOOST_REQUIRE(rv.first[s] != fwdpp::ts::TS_NULL_NODE);
        }
    BOOST_REQUIRE(tables.edges_are_minimally_sorted());
}

BOOST_AUTO_TEST_SUITE_END()
