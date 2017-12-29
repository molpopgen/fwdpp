/*!
  \file test_generalmut_vec.cc
  \ingroup unit
  \brief Testing fwdpp::generalmut_vec
*/
#include <unistd.h>
#include <config.h>
#include <sstream>
#include <boost/test/unit_test.hpp>
#include <fwdpp/diploid.hh>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/sugar/singlepop.hpp>
#include <fwdpp/sugar/infsites.hpp>
#include <fwdpp/sugar/generalmut.hpp>
#include <fwdpp/sugar/serialization.hpp>

struct generalmut_vec_tuple_wrapper
{
    fwdpp::generalmut_vec::constructor_tuple t;
    fwdpp::generalmut_vec m;
    static const std::vector<double> s, h;
    generalmut_vec_tuple_wrapper()
        : t(std::make_tuple(s, h, 2.2, 11, 17)), m(t)
    {
    }
};

const std::vector<double> generalmut_vec_tuple_wrapper::s = { { 0.0, -0.1 } };
const std::vector<double> generalmut_vec_tuple_wrapper::h = { { 1.0, 0.25 } };

BOOST_AUTO_TEST_SUITE(generalmut_vecTest)

BOOST_AUTO_TEST_CASE(construct_2)
{
    fwdpp::generalmut_vec p({ { 0.5, -1 } }, { { 1, 0 } }, 0.001, 1);
    BOOST_CHECK_EQUAL(p.s.size(), 2);
    BOOST_CHECK_EQUAL(p.h.size(), 2);
}

BOOST_AUTO_TEST_CASE(construct_2b)
{
    fwdpp::generalmut_vec p({ 0.5 }, { { 1, 0 } }, 0.001, 1);
    BOOST_CHECK_EQUAL(p.s.size(), 1);
    BOOST_CHECK_EQUAL(p.h.size(), 2);
}

BOOST_AUTO_TEST_CASE(construct_2c)
{
    fwdpp::generalmut_vec p({ 0.5 }, { 1 }, 0.001, 1);
    BOOST_CHECK_EQUAL(p.s.size(), 1);
    BOOST_CHECK_EQUAL(p.h.size(), 1);
}

BOOST_AUTO_TEST_CASE(construct_4)
{
    fwdpp::generalmut_vec p({ { 0.5, -1, 2.0, 3.0 } }, { { 1, 0, -1, 1 } },
                            0.001, 1);
    BOOST_CHECK_EQUAL(p.s.size(), 4);
    BOOST_CHECK_EQUAL(p.h.size(), 4);
}

BOOST_AUTO_TEST_CASE(construct_4b)
{
    fwdpp::generalmut_vec p({ 0.5 }, { { 1, 0, -1, 1 } }, 0.001, 1);
    BOOST_CHECK_EQUAL(p.s.size(), 1);
    BOOST_CHECK_EQUAL(p.h.size(), 4);
}

BOOST_AUTO_TEST_CASE(construct_4c)
{
    fwdpp::generalmut_vec p({ 0.5 }, { 1 }, 0.001, 1);
    BOOST_CHECK_EQUAL(p.s.size(), 1);
    BOOST_CHECK_EQUAL(p.h.size(), 1);
}

// Not implemented in library yet
BOOST_AUTO_TEST_CASE(serialize)
{
    fwdpp::generalmut_vec p({ { 0.5, -1 } }, { { 1, 0 } }, 0.001, 1);

    std::ostringstream o;
    fwdpp::io::serialize_mutation<fwdpp::generalmut_vec>()(p,o);

    std::istringstream i(o.str());
    auto p2 = fwdpp::io::deserialize_mutation<fwdpp::generalmut_vec>()(i);

    BOOST_CHECK_EQUAL(p.s.size(), p2.s.size());
    BOOST_CHECK_EQUAL(p.h.size(), p2.h.size());
    BOOST_CHECK_EQUAL(p.g, p2.g);
    BOOST_CHECK_EQUAL(p.pos, p2.pos);
}

BOOST_AUTO_TEST_CASE(copy_pop1)
{
    using mtype = fwdpp::generalmut_vec;
    using singlepop_t = fwdpp::singlepop<mtype>;
    singlepop_t pop1(100);
    singlepop_t pop2(pop1);
    BOOST_REQUIRE_EQUAL(pop1 == pop2, true);
}

BOOST_FIXTURE_TEST_CASE(construct_from_tuple, generalmut_vec_tuple_wrapper)
{
    BOOST_REQUIRE(m.s == s);
    BOOST_REQUIRE(m.h == h);
    BOOST_REQUIRE_EQUAL(m.pos, 2.2);
    BOOST_REQUIRE_EQUAL(m.g, 11);
    BOOST_REQUIRE_EQUAL(m.xtra, 17);
    BOOST_REQUIRE_EQUAL(m.neutral, false);
}

BOOST_AUTO_TEST_SUITE_END()
