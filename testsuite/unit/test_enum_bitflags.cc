#include <type_traits>
#include <boost/test/unit_test.hpp>
#include <fwdpp/util/enum_bitflags.hpp>

enum class E : unsigned
{
    a = 1 << 0,
    b = 1 << 1
};

BOOST_AUTO_TEST_SUITE(test_enum_bitflags)

BOOST_AUTO_TEST_CASE(test_zero_initialization)
{
    E e{};
    BOOST_REQUIRE_EQUAL(static_cast<std::underlying_type_t<E>>(e), 0);
}

BOOST_AUTO_TEST_CASE(test_bitwise_OR)
{
    using fwdpp::enum_bitflags::operator|;
    E e{E::a};
    BOOST_REQUIRE_EQUAL(!!(e | E::a), true);
    BOOST_REQUIRE_EQUAL(!!(e | E::b), true);
}

BOOST_AUTO_TEST_CASE(test_bitwise_OR_EQ)
{
    using fwdpp::enum_bitflags::operator|=;
    E e{E::a};
    e |= E::b;
    unsigned x = 1 << 0;
    x |= 1 << 1;
    BOOST_REQUIRE_EQUAL(static_cast<std::underlying_type_t<E>>(e), x);
}

BOOST_AUTO_TEST_CASE(test_bitwise_XOR)
{
    using fwdpp::enum_bitflags::operator^;
    E e{E::a};
    unsigned x = 1 << 0;
    x ^= 1 << 1;
    BOOST_REQUIRE_EQUAL(e^E::b, x);
}

BOOST_AUTO_TEST_CASE(test_bitwise_XOR_EQ)
{
    using fwdpp::enum_bitflags::operator^=;
    E e{E::a};
    e ^= E::b;
    unsigned x = 1 << 0;
    x ^= 1 << 1;
    BOOST_REQUIRE_EQUAL(static_cast<std::underlying_type_t<E>>(e), x);
}

BOOST_AUTO_TEST_CASE(test_bitwise_AND)
{
    using fwdpp::enum_bitflags::operator&;
    E e{E::a};
    BOOST_REQUIRE_EQUAL(!!(e & E::a), true);
    BOOST_REQUIRE_EQUAL(!!(e & E::b), false);
}

BOOST_AUTO_TEST_CASE(test_bitwise_AND_EQ)
{
    using fwdpp::enum_bitflags::operator&=;
    E e{E::a};
    e &= E::b;
    unsigned x = 1 << 0;
    x &= 1 << 1;
    BOOST_REQUIRE_EQUAL(static_cast<std::underlying_type_t<E>>(e), x);
}

BOOST_AUTO_TEST_CASE(test_bitwise_NEGATE)
{
    using fwdpp::enum_bitflags::operator~;
    E e{E::b};
    BOOST_REQUIRE_EQUAL(~e, ~2u);
}

BOOST_AUTO_TEST_SUITE_END()

