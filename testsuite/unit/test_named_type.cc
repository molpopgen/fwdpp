#include <vector>
#include <boost/test/unit_test.hpp>
#include <fwdpp/util/named_type.hpp>

BOOST_AUTO_TEST_SUITE(test_named_type)

BOOST_AUTO_TEST_CASE(test_move_construction)
{
    struct nt_tag{};
    using vec = fwdpp::strong_types::named_type<std::vector<int>, nt_tag>;
    std::vector<int> x{1,2,3,4,5,6};
    vec v{std::move(x)};
    BOOST_CHECK(x.empty());
}

BOOST_AUTO_TEST_SUITE_END()
