/*!
  \file type_traits_test.cc
  \ingroup unit
  \brief Testing fwdpp/type_traits.hpp

  These tests make sure that the type traits 
  actually return what we expect them to.
*/

#define BOOST_TEST_MODULE type_traits_test
#define BOOST_TEST_DYN_LINK

#include <config.h>
#include <boost/test/unit_test.hpp>
#include <fwdpp/diploid.hh>
#include <fwdpp/sugar/sugar.hpp>
using mutation_t = KTfwd::popgenmut;
using singlepop_t = KTfwd::singlepop<mutation_t>;

BOOST_AUTO_TEST_CASE( is_mmodel_t )
{
  singlepop_t pop(100);
  KTfwd::GSLrng_t<KTfwd::GSL_RNG_MT19937> r(101);
  auto mp = std::bind(KTfwd::infsites(),std::placeholders::_1,std::placeholders::_2,r.get(),std::ref(pop.mut_lookup),0,
		      0.001,0.,[&r](){return gsl_rng_uniform(r.get());},[](){return 0.;},[](){return 0.;});
  auto v = std::is_convertible<decltype(mp),KTfwd::traits::mmodel_t<singlepop_t::mcont_t> >::value;
  BOOST_REQUIRE_EQUAL(v,true);
}
