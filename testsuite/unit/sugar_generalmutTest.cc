/*!
  \file sugar_generalmutTest.cc
  \ingroup unit
  \brief Testing fwdpp::generalmut
*/
#include <unistd.h>
#include <config.h>
#include <sstream>
#include <boost/test/unit_test.hpp>
#include <fwdpp/sugar/singlepop.hpp>
#include <fwdpp/sugar/infsites.hpp>
#include <fwdpp/sugar/generalmut.hpp>

struct generalmut_tuple_wrapper
{
    fwdpp::generalmut<2>::constructor_tuple t;
    fwdpp::generalmut<2> m;
    static const std::array<double, 2> s, h;
    generalmut_tuple_wrapper() : t(std::make_tuple(s, h, 2.2, 11, 17)), m(t) {}
};

const std::array<double, 2> generalmut_tuple_wrapper::s = { { 0.0, -0.1 } };
const std::array<double, 2> generalmut_tuple_wrapper::h = { { 1.0, 0.25 } };

namespace fwdpp
{
    namespace io
    {
        template <> struct serialize_mutation<generalmut<4>>
        /// \brief Example specializing fwdpp::generalmut
        ///
        /// Every different size parameter is a new type,
        /// requiring a new specialization.  This example
        /// is for a size of 4.  The macro call below
        /// generates the call operator.
        {
            SPECIALIZE_SERIALIZE_MUTATION_GENERALMUT_BODY(4);
        };

        template <> struct deserialize_mutation<generalmut<4>>
        /// \brief Example specializing fwdpp::generalmut
        ///
        /// Every different size parameter is a new type,
        /// requiring a new specialization.  This example
        /// is for a size of 4.  The macro call below
        /// generates the call operator.
        {
            SPECIALIZE_DESERIALIZE_MUTATION_GENERALMUT_BODY(4);
        };
    }
}

BOOST_AUTO_TEST_SUITE(generalmutTest)

BOOST_AUTO_TEST_CASE(construct_2)
{
    fwdpp::generalmut<2> p({ { 0.5, -1 } }, { { 1, 0 } }, 0.001, 1);
    BOOST_CHECK_EQUAL(std::tuple_size<decltype(p)::array_t>(), 2);
    BOOST_CHECK_EQUAL(p.s.size(), 2);
    BOOST_CHECK_EQUAL(p.h.size(), 2);
}

BOOST_AUTO_TEST_CASE(construct_2b)
{
    fwdpp::generalmut<2> p({ { 0.5 } }, { { 1, 0 } }, 0.001, 1);
    BOOST_CHECK_EQUAL(std::tuple_size<decltype(p)::array_t>(), 2);
    BOOST_CHECK_EQUAL(p.s.size(), 2);
    BOOST_CHECK_EQUAL(p.h.size(), 2);
}

BOOST_AUTO_TEST_CASE(construct_2c)
{
    fwdpp::generalmut<2> p({ { 0.5 } }, { { 1 } }, 0.001, 1);
    BOOST_CHECK_EQUAL(std::tuple_size<decltype(p)::array_t>(), 2);
    BOOST_CHECK_EQUAL(p.s.size(), 2);
    BOOST_CHECK_EQUAL(p.h.size(), 2);
}

BOOST_AUTO_TEST_CASE(construct_4)
{
    fwdpp::generalmut<4> p({ { 0.5, -1, 2.0, 3.0 } }, { { 1, 0, -1, 1 } },
                           0.001, 1);
    BOOST_CHECK_EQUAL(std::tuple_size<decltype(p)::array_t>(), 4);
    BOOST_CHECK_EQUAL(p.s.size(), 4);
    BOOST_CHECK_EQUAL(p.h.size(), 4);
}

BOOST_AUTO_TEST_CASE(construct_4b)
{
    fwdpp::generalmut<4> p({ { 0.5 } }, { { 1, 0, -1, 1 } }, 0.001, 1);
    BOOST_CHECK_EQUAL(std::tuple_size<decltype(p)::array_t>(), 4);
    BOOST_CHECK_EQUAL(p.s.size(), 4);
    BOOST_CHECK_EQUAL(p.h.size(), 4);
}

BOOST_AUTO_TEST_CASE(construct_4c)
{
    fwdpp::generalmut<4> p({ { 0.5 } }, { { 1 } }, 0.001, 1);
    BOOST_CHECK_EQUAL(std::tuple_size<decltype(p)::array_t>(), 4);
    BOOST_CHECK_EQUAL(p.s.size(), 4);
    BOOST_CHECK_EQUAL(p.h.size(), 4);
}

BOOST_AUTO_TEST_CASE(serialize)
{
    fwdpp::generalmut<4> p({ { 0.5 } }, { { 1, 0, -1, 1 } }, 0.001, 1, 14);
    std::ostringstream o;
    fwdpp::io::serialize_mutation<fwdpp::generalmut<4>>()(p, o);
    std::istringstream i(o.str());
    auto p2 = fwdpp::io::deserialize_mutation<fwdpp::generalmut<4>>()(i);
    BOOST_REQUIRE(p == p2);
}

BOOST_AUTO_TEST_CASE(copy_pop1)
{
    using mtype = fwdpp::generalmut<2>;
    using singlepop_t = fwdpp::singlepop<mtype>;
    singlepop_t pop1(100);
    singlepop_t pop2(pop1);
    BOOST_REQUIRE_EQUAL(pop1 == pop2, true);
}

BOOST_FIXTURE_TEST_CASE(construct_from_tuple, generalmut_tuple_wrapper)
{
    BOOST_REQUIRE(m.s == s);
    BOOST_REQUIRE(m.h == h);
    BOOST_REQUIRE_EQUAL(m.pos, 2.2);
    BOOST_REQUIRE_EQUAL(m.g, 11);
    BOOST_REQUIRE_EQUAL(m.xtra, 17);
    BOOST_REQUIRE_EQUAL(m.neutral, false);
}

BOOST_AUTO_TEST_SUITE_END()
