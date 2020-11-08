#include <cmath>
#include <type_traits>
#include <boost/test/unit_test.hpp>
#include <fwdpp/util/validators.hpp>

BOOST_AUTO_TEST_SUITE(test_validators)

BOOST_AUTO_TEST_CASE(test_casts_to_int)
{
    BOOST_REQUIRE_NO_THROW({
        double x = 10.;
        fwdpp::validators::casts_to_int(x, "bad value");
    });

    BOOST_REQUIRE_THROW(
        {
            double x = 10.0000001;
            fwdpp::validators::casts_to_int(x, "bad value");
        },
        std::invalid_argument);

    BOOST_REQUIRE_THROW(
        {
            double x = std::nextafter(10., std::numeric_limits<double>::infinity());
            fwdpp::validators::casts_to_int(x, "bad value");
        },
        std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(test_non_negative)
{
    BOOST_REQUIRE_NO_THROW({
        double x = 10.;
        fwdpp::validators::non_negative(x, "bad value");
    });

    BOOST_REQUIRE_NO_THROW({
        double x = 0.;
        fwdpp::validators::non_negative(x, "bad value");
    });

    BOOST_REQUIRE_THROW(
        {
            double x = -10.;
            fwdpp::validators::non_negative(x, "bad value");
        },
        std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(test_is_positive)
{
    BOOST_REQUIRE_NO_THROW({
        double x = 10.;
        fwdpp::validators::is_positive(x, "bad value");
    });

    BOOST_REQUIRE_THROW(
        {
            double x = 0.;
            fwdpp::validators::is_positive(x, "bad value");
        },
        std::invalid_argument);

    BOOST_REQUIRE_THROW(
        {
            double x = -10.;
            fwdpp::validators::is_positive(x, "bad value");
        },
        std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(test_isfinite)
{
    BOOST_REQUIRE_NO_THROW({
        double x = 10.;
        fwdpp::validators::isfinite(x, "bad value");
    });

    BOOST_REQUIRE_THROW(
        {
            double x = 0. / 0.;
            fwdpp::validators::isfinite(x, "bad value");
        },
        std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()

