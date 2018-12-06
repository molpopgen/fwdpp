// This is really an API check

#include <config.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <fwdpp/wrapped_range.hpp>
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(test_wrapped_range)

BOOST_AUTO_TEST_CASE(test_summation)
{
    std::vector<int> x({ 1, 2, 3, 4, 5 });
    auto wr = fwdpp::make_wrapped_range(x.begin(), x.end());

    auto sum1 = std::accumulate(x.begin(), x.end(), 0);
    auto sum2 = std::accumulate(begin(wr), end(wr), 0);
    BOOST_REQUIRE_EQUAL(sum1, sum2);
    auto sum3 = 0;
    for (auto v : wr)
        {
            sum3 += v;
        }
    BOOST_REQUIRE_EQUAL(sum3, sum2);
}

BOOST_AUTO_TEST_SUITE_END()

