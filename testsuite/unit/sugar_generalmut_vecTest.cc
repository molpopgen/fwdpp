/*!
  \file sugar_generalmut_vecTest.cc
  \ingroup unit
  \brief Testing fwdpp::generalmut_vec
*/
#include <config.h>
#include <sstream>
#include <boost/test/unit_test.hpp>
#include <fwdpp/sugar/slocuspop.hpp>
#include <fwdpp/sugar/generalmut.hpp>

struct generalmut_vec_tuple_wrapper
{
    fwdpp::generalmut_vec::constructor_tuple t;
    fwdpp::generalmut_vec m;
    static const std::vector<std::tuple<double,double>> sh;
    generalmut_vec_tuple_wrapper()
        : t(std::make_tuple(sh, 2.2, 11, 17)), m(t)
    {
    }
};

using ttype = fwdpp::generalmut_vec::array_t::value_type;

const fwdpp::generalmut_vec::array_t generalmut_vec_tuple_wrapper::sh = { ttype{ 0.0,1.},ttype{ -0.1,0.25 } };

BOOST_AUTO_TEST_SUITE(generalmut_vecTest)

BOOST_AUTO_TEST_CASE(construct_2)
{
    fwdpp::generalmut_vec p({ ttype{ 0.5, -1 } , ttype { 1, 0 } }, 0.001, 1);
    BOOST_CHECK_EQUAL(p.sh.size(), 2);
}

BOOST_AUTO_TEST_CASE(serialize)
{
    fwdpp::generalmut_vec p({ ttype{ 0.5, -1 } ,  ttype{ 1, 0 } }, 0.001, 1);

    std::ostringstream o;
    fwdpp::io::serialize_mutation<fwdpp::generalmut_vec>()(o, p);

    std::istringstream i(o.str());
    auto p2 = fwdpp::io::deserialize_mutation<fwdpp::generalmut_vec>()(i);

    BOOST_REQUIRE(p.sh == p2.sh);
    BOOST_CHECK_EQUAL(p.sh.size(), p2.sh.size());
    BOOST_CHECK_EQUAL(p.g, p2.g);
    BOOST_CHECK_EQUAL(p.pos, p2.pos);
}

BOOST_AUTO_TEST_CASE(copy_pop1)
{
    using mtype = fwdpp::generalmut_vec;
    using slocuspop_t = fwdpp::slocuspop<mtype>;
    slocuspop_t pop1(100);
    slocuspop_t pop2(pop1);
    BOOST_REQUIRE_EQUAL(pop1 == pop2, true);
}

BOOST_FIXTURE_TEST_CASE(construct_from_tuple, generalmut_vec_tuple_wrapper)
{
    BOOST_REQUIRE(m.sh == sh);
    BOOST_REQUIRE_EQUAL(m.pos, 2.2);
    BOOST_REQUIRE_EQUAL(m.g, 11);
    BOOST_REQUIRE_EQUAL(m.xtra, 17);
    BOOST_REQUIRE_EQUAL(m.neutral, false);
}

BOOST_AUTO_TEST_SUITE_END()
