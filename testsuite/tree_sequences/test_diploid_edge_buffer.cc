#include <boost/test/unit_test.hpp>
#include <gsl/gsl_randist.h>
#include <fwdpp/GSLrng_t.hpp>
#include <fwdpp/ts/diploid_edge_buffer.hpp>

std::pair<fwdpp::ts::TS_NODE_INT, unsigned>
diploid_wf(const fwdpp::GSLrng_mt& rng, unsigned N,
           unsigned parental_generation, unsigned ngenerations,
           fwdpp::ts::TS_NODE_INT first_parental_index,
           fwdpp::ts::table_collection& tables,
           fwdpp::ts::diploid_edge_buffer& buffer)
// Based on the tskit tutorial examples https://tskit-dev.github.io/tutorials/wfforward.html
{
    decltype(first_parental_index) parent_index_offset = first_parental_index;
    fwdpp::ts::TS_NODE_INT next_offspring_node;
    std::vector<double> breakpoints; //No rec.
    const fwdpp::ts::TS_NODE_INT offspring_deme{ 0 };
    unsigned generation = 1;
    for (; generation <= ngenerations; ++generation)
        {
            buffer.begin_epoch(N);
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
                    next_offspring_node
                        = fwdpp::ts::register_offspring_and_buffer_edges(
                            breakpoints, parent, std::make_tuple(pn0, pn1),
                            offspring_deme, parental_generation + generation,
                            tables, buffer);
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

                    next_offspring_node
                        = fwdpp::ts::register_offspring_and_buffer_edges(
                            breakpoints, parent, std::make_tuple(pn0, pn1),
                            offspring_deme, parental_generation + generation,
                            tables, buffer);
                    if (next_offspring_node == fwdpp::ts::TS_NULL_NODE)
                        {
                            throw std::runtime_error("NULL offspring node");
                        }
                }
            parent_index_offset += 2 * N;
            buffer.end_epoch();
        }
    return std::make_pair(parent_index_offset,
                          parental_generation + generation - 1);
}

struct wfnorec
{
    const double L;
    const fwdpp::ts::TS_NODE_INT N;
    fwdpp::ts::table_collection tables;
    fwdpp::ts::diploid_edge_buffer buffer;
    fwdpp::GSLrng_mt rng;
    unsigned parental_generation;
    wfnorec()
        : L{ 1.0 }, N{ 1000 },
          tables(2 * N, 0.0, 0, L), buffer{}, rng{ 42 }, parental_generation{
              0
          }
    {
    }
};

BOOST_AUTO_TEST_SUITE(test_diploid_edge_buffer)

BOOST_FIXTURE_TEST_CASE(test_sort_order, wfnorec)
{
    auto p = diploid_wf(rng, N, parental_generation, 100, 0, tables, buffer);
    parental_generation = p.second;
    buffer.prepare_tables_for_simplification(tables);
    BOOST_REQUIRE_EQUAL(tables.edge_table.size(), parental_generation * 2 * N);
    BOOST_REQUIRE_EQUAL(tables.node_table.size(),
                        (parental_generation + 1) * 2 * N);
    BOOST_REQUIRE_EQUAL(tables.edges_are_sorted(), true);
}

BOOST_FIXTURE_TEST_CASE(test_two_rounds, wfnorec)
{
    auto p = diploid_wf(rng, N, parental_generation, 100, 0, tables, buffer);
    buffer.prepare_tables_for_simplification(tables);
    tables.edge_offset = tables.edge_table.size();
    parental_generation = p.second;
    p = diploid_wf(rng, N, parental_generation, 100, p.first, tables, buffer);
    parental_generation = p.second;
    buffer.prepare_tables_for_simplification(tables);
    BOOST_REQUIRE_EQUAL(tables.edge_table.size(), parental_generation * 2 * N);
    BOOST_REQUIRE_EQUAL(tables.node_table.size(),
                        (parental_generation + 1) * 2 * N);
    BOOST_REQUIRE_EQUAL(tables.edges_are_sorted(), true);
}

BOOST_AUTO_TEST_SUITE_END()
