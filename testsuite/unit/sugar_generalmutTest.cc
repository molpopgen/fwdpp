/*!
  \file sugar_generalmutTest.cc
  \ingroup unit
  \brief Testing KTfwd::generalmut
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

struct generalmut_tuple_wrapper
{
    KTfwd::generalmut<2>::constructor_tuple t;
    KTfwd::generalmut<2> m;
    static const std::array<double, 2> s, h;
    generalmut_tuple_wrapper() : t(std::make_tuple(s, h, 2.2, 11, 17)), m(t) {}
};

const std::array<double, 2> generalmut_tuple_wrapper::s = { { 0.0, -0.1 } };
const std::array<double, 2> generalmut_tuple_wrapper::h = { { 1.0, 0.25 } };

BOOST_AUTO_TEST_SUITE(generalmutTest)

BOOST_AUTO_TEST_CASE(construct_2)
{
    KTfwd::generalmut<2> p({ { 0.5, -1 } }, { { 1, 0 } }, 0.001, 1);
    BOOST_CHECK_EQUAL(std::tuple_size<decltype(p)::array_t>(), 2);
    BOOST_CHECK_EQUAL(p.s.size(), 2);
    BOOST_CHECK_EQUAL(p.h.size(), 2);
}

BOOST_AUTO_TEST_CASE(construct_2b)
{
    KTfwd::generalmut<2> p({ { 0.5 } }, { { 1, 0 } }, 0.001, 1);
    BOOST_CHECK_EQUAL(std::tuple_size<decltype(p)::array_t>(), 2);
    BOOST_CHECK_EQUAL(p.s.size(), 2);
    BOOST_CHECK_EQUAL(p.h.size(), 2);
}

BOOST_AUTO_TEST_CASE(construct_2c)
{
    KTfwd::generalmut<2> p({ { 0.5 } }, { { 1 } }, 0.001, 1);
    BOOST_CHECK_EQUAL(std::tuple_size<decltype(p)::array_t>(), 2);
    BOOST_CHECK_EQUAL(p.s.size(), 2);
    BOOST_CHECK_EQUAL(p.h.size(), 2);
}

BOOST_AUTO_TEST_CASE(construct_4)
{
    KTfwd::generalmut<4> p({ { 0.5, -1, 2.0, 3.0 } }, { { 1, 0, -1, 1 } },
                           0.001, 1);
    BOOST_CHECK_EQUAL(std::tuple_size<decltype(p)::array_t>(), 4);
    BOOST_CHECK_EQUAL(p.s.size(), 4);
    BOOST_CHECK_EQUAL(p.h.size(), 4);
}

BOOST_AUTO_TEST_CASE(construct_4b)
{
    KTfwd::generalmut<4> p({ { 0.5 } }, { { 1, 0, -1, 1 } }, 0.001, 1);
    BOOST_CHECK_EQUAL(std::tuple_size<decltype(p)::array_t>(), 4);
    BOOST_CHECK_EQUAL(p.s.size(), 4);
    BOOST_CHECK_EQUAL(p.h.size(), 4);
}

BOOST_AUTO_TEST_CASE(construct_4c)
{
    KTfwd::generalmut<4> p({ { 0.5 } }, { { 1 } }, 0.001, 1);
    BOOST_CHECK_EQUAL(std::tuple_size<decltype(p)::array_t>(), 4);
    BOOST_CHECK_EQUAL(p.s.size(), 4);
    BOOST_CHECK_EQUAL(p.h.size(), 4);
}

BOOST_AUTO_TEST_CASE(serialize)
{
    KTfwd::generalmut<2> p({ { 0.5, -1 } }, { { 1, 0 } }, 0.001, 1);

    std::ostringstream o;
    KTfwd::mutation_writer w;
    w(p, o);

    KTfwd::mutation_reader<decltype(p)> r;
    std::istringstream i(o.str());
    auto p2 = r(i);

    BOOST_CHECK_EQUAL(p.s.size(), p2.s.size());
    BOOST_CHECK_EQUAL(p.h.size(), p2.h.size());
    BOOST_CHECK_EQUAL(p.g, p2.g);
    BOOST_CHECK_EQUAL(p.pos, p2.pos);
}

BOOST_AUTO_TEST_CASE(serialize_gz)
{
    KTfwd::generalmut<2> p({ { 0.5, -1 } }, { { 1, 0 } }, 0.001, 1);

    gzFile out = gzopen("test_generalmut_file.gz", "w");
    KTfwd::mutation_writer w;
    w(p, out);
    gzclose(out);

    KTfwd::mutation_reader<decltype(p)> r;
    out = gzopen("test_generalmut_file.gz", "r");
    auto p2 = r(out);

    BOOST_CHECK_EQUAL(p.s.size(), p2.s.size());
    BOOST_CHECK_EQUAL(p.h.size(), p2.h.size());
    BOOST_CHECK_EQUAL(p.g, p2.g);
    BOOST_CHECK_EQUAL(p.pos, p2.pos);

    unlink("test_generalmut_file.gz");
}

BOOST_AUTO_TEST_CASE(serialize_pop1)
{
    using mtype = KTfwd::generalmut<2>;
    using singlepop_t = KTfwd::singlepop<mtype>;
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
