#include <config.h>
#include <limits>
#include <cstdint>
#include <iostream>
#include <fwdpp/util/nested_forward_lists.hpp>
#include <boost/test/unit_test.hpp>

namespace
{
    struct Foo
    {
    };
}

BOOST_AUTO_TEST_SUITE(test_nested_forward_lists)

BOOST_AUTO_TEST_CASE(test_signed_overflow)
{
    using buffer_t = fwdpp::nested_forward_lists<Foo, std::int8_t, -1>;
    buffer_t buffer;
    buffer.reset(1);
    auto m = std::numeric_limits<std::int8_t>::max();
    BOOST_REQUIRE_THROW(
        {
            for (decltype(m) i = 0; i <= m; ++i)
                {
                    buffer.extend(0, Foo{});
                }
        },
        fwdpp::nested_forward_lists_overflow);
}

BOOST_AUTO_TEST_CASE(test_unsigned_overflow)
{
    using buffer_t
        = fwdpp::nested_forward_lists<Foo, std::uint8_t,
                                      std::numeric_limits<std::uint8_t>::max()>;
    buffer_t buffer;
    buffer.reset(1);
    auto m = std::numeric_limits<std::int8_t>::max();
    BOOST_REQUIRE_THROW(
        {
            for (decltype(m) i = 0; i <= m; ++i)
                {
                    buffer.extend(0, Foo{});
                }
        },
        fwdpp::nested_forward_lists_overflow);
}

BOOST_AUTO_TEST_SUITE_END()

