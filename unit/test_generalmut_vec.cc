/*! 
  \file test_generalmut_vec.cc 
  \ingroup unit 
  \brief Testing KTfwd::generalmut_vec
*/
#define BOOST_TEST_MODULE generalmut_vecTest
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
  KTfwd::generalmut_vec p( {{0.5,-1}},{{1,0}},0.001,1,2);
  BOOST_CHECK_EQUAL(p.s.size(),2);
  BOOST_CHECK_EQUAL(p.h.size(),2);
}

BOOST_AUTO_TEST_CASE( construct_2b )
{
  KTfwd::generalmut_vec p( {0.5},{{1,0}},0.001,1,2);
  BOOST_CHECK_EQUAL(p.s.size(),1);
  BOOST_CHECK_EQUAL(p.h.size(),2);
}

BOOST_AUTO_TEST_CASE( construct_2c )
{
  KTfwd::generalmut_vec p( {0.5},{1},0.001,1,2);
  BOOST_CHECK_EQUAL(p.s.size(),1);
  BOOST_CHECK_EQUAL(p.h.size(),1);
}

BOOST_AUTO_TEST_CASE( construct_4 )
{
  KTfwd::generalmut_vec p( {{0.5,-1,2.0,3.0}},{{1,0,-1,1}},0.001,1,2);
  BOOST_CHECK_EQUAL(p.s.size(),4);
  BOOST_CHECK_EQUAL(p.h.size(),4);
}

BOOST_AUTO_TEST_CASE( construct_4b )
{
  KTfwd::generalmut_vec p( {0.5},{{1,0,-1,1}},0.001,1,2);
  BOOST_CHECK_EQUAL(p.s.size(),1);
  BOOST_CHECK_EQUAL(p.h.size(),4);
}

BOOST_AUTO_TEST_CASE( construct_4c )
{
  KTfwd::generalmut_vec p( {0.5},{1},0.001,1,2);
  BOOST_CHECK_EQUAL(p.s.size(),1);
  BOOST_CHECK_EQUAL(p.h.size(),1);
}

//Not implemented in library yet
BOOST_AUTO_TEST_CASE( serialize )
{
  KTfwd::generalmut_vec p( {{0.5,-1}},{{1,0}},0.001,1,2);

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
  KTfwd::generalmut_vec p( {{0.5,-1}},{{1,0}},0.001,1,2);

  gzFile out = gzopen("test_generalmut_vec_file.gz","w");
  KTfwd::mutation_writer w;
  w(p,out);
  gzclose(out);
  
  KTfwd::mutation_reader<decltype(p)> r;
  out = gzopen("test_generalmut_vec_file.gz","r");
  auto p2 = r(out);

  BOOST_CHECK_EQUAL(p.s.size(),p2.s.size());
  BOOST_CHECK_EQUAL(p.h.size(),p2.h.size());
  BOOST_CHECK_EQUAL(p.n,p2.n);
  BOOST_CHECK_EQUAL(p.g,p2.g);
  BOOST_CHECK_EQUAL(p.pos,p2.pos);

  unlink("test_generalmut_vec_file.gz");
}
