#include <iostream>
#include <boost/test/unit_test.hpp>
#include <fwdpp/ts/simplification/simplification.hpp>

namespace
{
    std::int64_t
    get_list_length(const fwdpp::ts::simplification::ancestry_list& al, std::size_t i)
    {
        int len = 0;
        auto f = al.head(i);

        while (f != -1)
            {
                ++len;
                f = al.next(f);
            }
        return len;
    }
} // namespace

BOOST_AUTO_TEST_SUITE(test_ancestry_list)

BOOST_AUTO_TEST_CASE(test_construct_and_fill)
// The data added don't refer to valid
// tree sequences and are just for testing
{
    fwdpp::ts::simplification::ancestry_list al;
    al.reset(4);
    al.extend(0, 0., 1., 3);
    al.extend(1, 0., 0.5, 3);
    al.extend(0, 0., 0.5, 4);

    BOOST_REQUIRE_EQUAL(get_list_length(al, 0), 2);
    BOOST_REQUIRE_EQUAL(get_list_length(al, 1), 1);
}

BOOST_AUTO_TEST_CASE(test_fill_nullify_once_fill)
{
    fwdpp::ts::simplification::ancestry_list al;
    al.reset(4);
    al.extend(0, 0., 1., 3);
    al.extend(1, 0., 0.5, 3);
    al.extend(0, 0., 0.5, 4);
    BOOST_REQUIRE_EQUAL(get_list_length(al, 0), 2);
    BOOST_REQUIRE_EQUAL(al.tail(0), 2);
    al.nullify_list(0);
    BOOST_REQUIRE_EQUAL(get_list_length(al, 0), 0);
    BOOST_REQUIRE_EQUAL(al.tail(0), -1);

    al.extend(0, 0., 1., 2);
    BOOST_REQUIRE_EQUAL(get_list_length(al, 0), 1);
}

BOOST_AUTO_TEST_SUITE_END()

