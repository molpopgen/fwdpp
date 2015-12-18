/*! 
  \file test_generalmut.cc 
  \ingroup unit 
  \brief Testing KTfwd::generalmut
*/
#define BOOST_TEST_MODULE generalmutTest
#define BOOST_TEST_DYN_LINK 

#include <unistd.h>
#include <config.h>
#include <sstream>
#include <boost/test/unit_test.hpp>
#include <fwdpp/diploid.hh>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/sugar/singlepop.hpp>
#include <fwdpp/sugar/infsites.hpp>
#include <fwdpp/sugar/generalmut.hpp>

BOOST_AUTO_TEST_CASE( construct_2 )
{
  KTfwd::generalmut<2> p( {{0.5,-1}},{{1,0}},0.001,1,2);
  BOOST_CHECK_EQUAL(std::tuple_size<decltype(p)::array_t>(),2);
  BOOST_CHECK_EQUAL(p.s.size(),2);
  BOOST_CHECK_EQUAL(p.h.size(),2);
}

BOOST_AUTO_TEST_CASE( construct_2b )
{
  KTfwd::generalmut<2> p( {{0.5}},{{1,0}},0.001,1,2);
  BOOST_CHECK_EQUAL(std::tuple_size<decltype(p)::array_t>(),2);
  BOOST_CHECK_EQUAL(p.s.size(),2);
  BOOST_CHECK_EQUAL(p.h.size(),2);
}

BOOST_AUTO_TEST_CASE( construct_2c )
{
  KTfwd::generalmut<2> p( {{0.5}},{{1}},0.001,1,2);
  BOOST_CHECK_EQUAL(std::tuple_size<decltype(p)::array_t>(),2);
  BOOST_CHECK_EQUAL(p.s.size(),2);
  BOOST_CHECK_EQUAL(p.h.size(),2);
}

BOOST_AUTO_TEST_CASE( construct_4 )
{
  KTfwd::generalmut<4> p( {{0.5,-1,2.0,3.0}},{{1,0,-1,1}},0.001,1,2);
  BOOST_CHECK_EQUAL(std::tuple_size<decltype(p)::array_t>(),4);
  BOOST_CHECK_EQUAL(p.s.size(),4);
  BOOST_CHECK_EQUAL(p.h.size(),4);
}

BOOST_AUTO_TEST_CASE( construct_4b )
{
  KTfwd::generalmut<4> p( {{0.5}},{{1,0,-1,1}},0.001,1,2);
  BOOST_CHECK_EQUAL(std::tuple_size<decltype(p)::array_t>(),4);
  BOOST_CHECK_EQUAL(p.s.size(),4);
  BOOST_CHECK_EQUAL(p.h.size(),4);
}

BOOST_AUTO_TEST_CASE( construct_4c )
{
  KTfwd::generalmut<4> p( {{0.5}},{{1}},0.001,1,2);
  BOOST_CHECK_EQUAL(std::tuple_size<decltype(p)::array_t>(),4);
  BOOST_CHECK_EQUAL(p.s.size(),4);
  BOOST_CHECK_EQUAL(p.h.size(),4);
}

BOOST_AUTO_TEST_CASE( serialize )
{
  KTfwd::generalmut<2> p( {{0.5,-1}},{{1,0}},0.001,1,2);

  std::ostringstream o;
  KTfwd::mutation_writer w;
  w(p,o);

  KTfwd::mutation_reader<decltype(p)> r;
  std::istringstream i(o.str());
  auto p2 = r(i);

  BOOST_CHECK_EQUAL(p.s.size(),p2.s.size());
  BOOST_CHECK_EQUAL(p.h.size(),p2.h.size());
  BOOST_CHECK_EQUAL(p.n,p2.n);
  BOOST_CHECK_EQUAL(p.g,p2.g);
  BOOST_CHECK_EQUAL(p.pos,p2.pos);
}

BOOST_AUTO_TEST_CASE( serialize_gz )
{
  KTfwd::generalmut<2> p( {{0.5,-1}},{{1,0}},0.001,1,2);

  gzFile out = gzopen("test_generalmut_file.gz","w");
  KTfwd::mutation_writer w;
  w(p,out);
  gzclose(out);
  
  KTfwd::mutation_reader<decltype(p)> r;
  out = gzopen("test_generalmut_file.gz","r");
  auto p2 = r(out);

  BOOST_CHECK_EQUAL(p.s.size(),p2.s.size());
  BOOST_CHECK_EQUAL(p.h.size(),p2.h.size());
  BOOST_CHECK_EQUAL(p.n,p2.n);
  BOOST_CHECK_EQUAL(p.g,p2.g);
  BOOST_CHECK_EQUAL(p.pos,p2.pos);

  unlink("test_generalmut_file.gz");
}

BOOST_AUTO_TEST_CASE( serialize_pop1 )
{
  using mtype = KTfwd::generalmut<2>;
  using singlepop_serialized_t = KTfwd::singlepop_serialized<mtype,
							     KTfwd::mutation_writer,
							     KTfwd::mutation_reader<mtype>>;

  singlepop_serialized_t pop1(100);
  singlepop_serialized_t pop2(pop1);
  
}

BOOST_AUTO_TEST_CASE( construct_extinct )
{
  auto x = KTfwd::generalmut<2>(KTfwd::tags::extinct());
  BOOST_CHECK_EQUAL(x.n,0);
  BOOST_CHECK_EQUAL(x.g,std::numeric_limits<unsigned>::max());
}
