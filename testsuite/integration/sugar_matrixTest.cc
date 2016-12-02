/* \brief Integration and unit tests of data matrix generation
 * \ingroup unit
 */

#include <config.h>
#include <boost/test/unit_test.hpp>
#include "../fixtures/sugar_fixtures.hpp"
#include "../util/quick_evolve_sugar.hpp"
#include <fwdpp/sugar/matrix.hpp>
#include <fwdpp/sugar/sampling.hpp>
#include <gsl/gsl_matrix_short.h>

BOOST_AUTO_TEST_CASE(singlepop_hapmatrix)
{
    using spoptype = singlepop_popgenmut_fixture::poptype;
    spoptype pop(1000);
    simulate_singlepop(pop, 10000);
    auto keys = mutation_keys(pop, { 0, 1, 2, 3 }, true, true);
    auto m = haplotype_matrix(pop, { 0, 1, 2, 3 }, keys.first, keys.second);

    BOOST_REQUIRE_EQUAL(m.nrow, 8);
    BOOST_REQUIRE_EQUAL(m.neutral.size(), m.nrow * keys.first.size());
    BOOST_REQUIRE_EQUAL(m.selected.size(), m.nrow * keys.second.size());
    BOOST_REQUIRE_EQUAL(keys.first.size(), m.neutral_positions.size());
    BOOST_REQUIRE_EQUAL(keys.first.size(), m.neutral_popfreq.size());
    BOOST_REQUIRE_EQUAL(keys.second.size(), m.selected_positions.size());
    BOOST_REQUIRE_EQUAL(keys.second.size(), m.selected_popfreq.size());
}

BOOST_AUTO_TEST_CASE(singlepop_hapmatrix_compare_to_sample)
{
    using spoptype = singlepop_popgenmut_fixture::poptype;
    spoptype pop(1000);
    simulate_singlepop(pop, 10000);
    std::vector<std::size_t> indlist;
    // Sample a LOT of individuals
    for (std::size_t i = 100; i < 750; i += 5)
        indlist.push_back(i);
    std::vector<unsigned> indlist2(indlist.begin(), indlist.end());
    auto keys = mutation_keys(pop, indlist, true, true);
    // keys come out unsorted, so we have to sort for this test:
    std::sort(keys.first.begin(), keys.first.end(),
              [&pop](const std::pair<std::size_t, KTfwd::uint_t> &a,
                     const std::pair<std::size_t, KTfwd::uint_t> &b) {
                  return pop.mutations[a.first].pos
                         < pop.mutations[b.first].pos;
              });
    std::sort(keys.second.begin(), keys.second.end(),
              [&pop](const std::pair<std::size_t, KTfwd::uint_t> &a,
                     const std::pair<std::size_t, KTfwd::uint_t> &b) {
                  return pop.mutations[a.first].pos
                         < pop.mutations[b.first].pos;
              });
    auto m = haplotype_matrix(pop, indlist, keys.first, keys.second);
    // Same basic checks as previous tests
    BOOST_REQUIRE_EQUAL(m.nrow, 2 * indlist.size());
    BOOST_REQUIRE_EQUAL(m.neutral.size(), m.nrow * keys.first.size());
    BOOST_REQUIRE_EQUAL(m.selected.size(), m.nrow * keys.second.size());
    BOOST_REQUIRE_EQUAL(keys.first.size(), m.neutral_positions.size());
    BOOST_REQUIRE_EQUAL(keys.first.size(), m.neutral_popfreq.size());
    BOOST_REQUIRE_EQUAL(keys.second.size(), m.selected_positions.size());
    BOOST_REQUIRE_EQUAL(keys.second.size(), m.selected_popfreq.size());

    // Now, compare to an independent calculation of the genotypes
    // for same individuals
    auto s = KTfwd::sample_separate(pop, indlist2, true);
    // Check same total # of mutations
    BOOST_REQUIRE_EQUAL(s.first.size(), m.neutral_positions.size());
    BOOST_REQUIRE_EQUAL(s.second.size(), m.selected_positions.size());
    std::vector<double> pos;
    for (auto &&i : s.first)
        pos.push_back(i.first);
    // Make sure all positions come out ok
    for (std::size_t i = 0; i < pos.size(); ++i)
        {
            BOOST_REQUIRE_EQUAL(pos[i], m.neutral_positions[i]);
        }
    pos.clear();
    for (auto &&i : s.second)
        pos.push_back(i.first);
    for (std::size_t i = 0; i < pos.size(); ++i)
        {
            BOOST_REQUIRE_EQUAL(pos[i], m.selected_positions[i]);
        }
    // Use GSL matrix stuff to make sure that matrix == s
    auto v = gsl_matrix_short_const_view_array(m.neutral.data(), m.nrow,
                                               m.neutral_positions.size());
    BOOST_REQUIRE_EQUAL(v.matrix.size1, 2 * indlist.size());
    // Make sure column sums check out
    for (std::size_t c = 0; c < v.matrix.size2; ++c)
        {
            auto col_view = gsl_matrix_short_const_column(&v.matrix, c);
            unsigned sum = 0, sum2 = 0;
            for (std::size_t i = 0; i < v.matrix.size1; ++i)
                {
                    sum += gsl_vector_short_get(&col_view.vector, i);
                    sum2 += (s.first[c].second[i] == '1') ? 1 : 0;
                }
            BOOST_REQUIRE_EQUAL(sum, sum2);
        }
    // Now, row sums.
    for (std::size_t r = 0; r < v.matrix.size1; ++r)
        {
            auto row_view = gsl_matrix_short_const_row(&v.matrix, r);
            unsigned sum = 0, sum2 = 0;
            for (std::size_t i = 0; i < v.matrix.size2; ++i)
                {
                    sum += gsl_vector_short_get(&row_view.vector, i);
                    sum2 += (s.first[i].second[r] == '1') ? 1 : 0;
                }
            BOOST_REQUIRE_EQUAL(sum, sum2);
        }
    // Repeat for selected sites
    if (!m.selected_positions.empty())
        {
            v = gsl_matrix_short_const_view_array(m.selected.data(), m.nrow,
                                                  m.selected_positions.size());
            BOOST_REQUIRE_EQUAL(v.matrix.size1, 2 * indlist.size());
            // Make sure column sums check out
            for (std::size_t c = 0; c < v.matrix.size2; ++c)
                {
                    auto col_view
                        = gsl_matrix_short_const_column(&v.matrix, c);
                    unsigned sum = 0, sum2 = 0;
                    for (std::size_t i = 0; i < v.matrix.size1; ++i)
                        {
                            sum += gsl_vector_short_get(&col_view.vector, i);
                            sum2 += (s.second[c].second[i] == '1') ? 1 : 0;
                        }
                    BOOST_REQUIRE_EQUAL(sum, sum2);
                }
            // Now, row sums.
            for (std::size_t r = 0; r < v.matrix.size1; ++r)
                {
                    auto row_view = gsl_matrix_short_const_row(&v.matrix, r);
                    unsigned sum = 0, sum2 = 0;
                    for (std::size_t i = 0; i < v.matrix.size2; ++i)
                        {
                            sum += gsl_vector_short_get(&row_view.vector, i);
                            sum2 += (s.second[i].second[r] == '1') ? 1 : 0;
                        }
                    BOOST_REQUIRE_EQUAL(sum, sum2);
                }
        }
}

BOOST_AUTO_TEST_CASE(singlepop_hapmatrix_GSL_behavior)
// What happens if try to create empty matrix?
{
    using spoptype = singlepop_popgenmut_fixture::poptype;
    spoptype pop(1000);
    std::vector<std::size_t> indlist;
    // Sample a LOT of individuals
    for (std::size_t i = 100; i < 750; i += 5)
        indlist.push_back(i);
    auto keys = mutation_keys(pop, indlist, true, true);
    auto m = haplotype_matrix(pop, indlist, keys.first, keys.second);
    BOOST_REQUIRE(keys.first.empty());
    // Use GSL matrix stuff to make sure that matrix == s
    auto v = gsl_matrix_short_const_view_array(m.neutral.data(), m.nrow,
                                               m.neutral_positions.size());
}

BOOST_AUTO_TEST_CASE(singlepop_genotype_matrix)
{
    using spoptype = singlepop_popgenmut_fixture::poptype;
    spoptype pop(1000);
    simulate_singlepop(pop, 10000);
    auto keys = mutation_keys(pop, { 0, 1, 2, 3 }, true, true);
    auto m = genotype_matrix(pop, { 0, 1, 2, 3 }, keys.first, keys.second);

    BOOST_REQUIRE_EQUAL(m.nrow, 4);
    BOOST_REQUIRE_EQUAL(m.neutral.size(), m.nrow * keys.first.size());
    BOOST_REQUIRE_EQUAL(m.selected.size(), m.nrow * keys.second.size());
    BOOST_REQUIRE_EQUAL(keys.first.size(), m.neutral_positions.size());
    BOOST_REQUIRE_EQUAL(keys.first.size(), m.neutral_popfreq.size());
    BOOST_REQUIRE_EQUAL(keys.second.size(), m.selected_positions.size());
    BOOST_REQUIRE_EQUAL(keys.second.size(), m.selected_popfreq.size());
}
