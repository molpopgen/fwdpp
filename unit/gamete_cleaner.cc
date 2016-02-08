/*!
  \file gamete_cleaner.cc
  \ingroup unit
  \brief Testing KTfwd::fwdpp_internal::gamete_cleaner
*/
#define BOOST_TEST_MODULE sugar_singlepop
#define BOOST_TEST_DYN_LINK

#include <config.h>
#include <boost/test/unit_test.hpp>
#include <fwdpp/diploid.hh>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/sugar/singlepop.hpp>
#include <fwdpp/sugar/infsites.hpp>

using mutation_t = KTfwd::popgenmut;
using singlepop_t = KTfwd::singlepop<mutation_t>;

BOOST_AUTO_TEST_CASE( test_remove_all )
{
  singlepop_t pop(1000);
  pop.mutations.emplace_back( mutation_t(0,0,0,0) );
  pop.mutations.emplace_back( mutation_t(1,0,0,0) );
  pop.mcounts.emplace_back(2000); //fixed
  pop.mcounts.emplace_back(1); //singleton
  pop.gametes[0].mutations.push_back(0);
  pop.gametes[0].n=1999;
  pop.gametes.push_back( singlepop_t::gamete_t(1) );
  pop.gametes[1].mutations.push_back(1);
  BOOST_REQUIRE_EQUAL( KTfwd::check_sum(pop.gametes,2000), true );
  BOOST_REQUIRE_EQUAL( KTfwd::popdata_sane(pop.diploids,pop.gametes,pop.mcounts), true );
  KTfwd::fwdpp_internal::gamete_cleaner(pop.gametes,pop.mutations,pop.mcounts,2000,std::true_type());
  BOOST_REQUIRE_EQUAL( pop.gametes[0].mutations.size(),0 );
}

BOOST_AUTO_TEST_CASE( test_remove_nothing )
{
  singlepop_t pop(1000);
  pop.mutations.emplace_back( mutation_t(0,0,0,0) );
  pop.mutations.emplace_back( mutation_t(1,0,0,0) );
  pop.mcounts.emplace_back(2000); //fixed
  pop.mcounts.emplace_back(1); //singleton
  pop.gametes[0].mutations.push_back(0);
  pop.gametes[0].n=1999;
  pop.gametes.push_back( singlepop_t::gamete_t(1) );
  pop.gametes[1].mutations.push_back(1);
  BOOST_REQUIRE_EQUAL( KTfwd::check_sum(pop.gametes,2000), true );
  BOOST_REQUIRE_EQUAL( KTfwd::popdata_sane(pop.diploids,pop.gametes,pop.mcounts), true );
  KTfwd::fwdpp_internal::gamete_cleaner(pop.gametes,pop.mutations,pop.mcounts,2000,KTfwd::remove_nothing());
  BOOST_REQUIRE_EQUAL( pop.gametes[0].mutations.size(),1 );
}

BOOST_AUTO_TEST_CASE( test_remove_neutral )
{
  singlepop_t pop(1000);
  pop.mutations.emplace_back( mutation_t(0,0,0,0) );
  pop.mutations.emplace_back( mutation_t(1,1,0,0) ); //not neutral
  pop.mcounts.emplace_back(2000); //fixed
  pop.mcounts.emplace_back(2000); //singleton
  pop.gametes[0].mutations.push_back(0);
  pop.gametes[0].mutations.push_back(1);
  BOOST_REQUIRE_EQUAL( KTfwd::check_sum(pop.gametes,2000), true );
  BOOST_REQUIRE_EQUAL( KTfwd::popdata_sane(pop.diploids,pop.gametes,pop.mcounts), true );
  KTfwd::fwdpp_internal::gamete_cleaner(pop.gametes,pop.mutations,pop.mcounts,2000,KTfwd::remove_neutral());
  BOOST_REQUIRE_EQUAL( pop.gametes[0].mutations.size(),1 );
  BOOST_REQUIRE_EQUAL( pop.mutations[pop.gametes[0].mutations[0]].neutral,false );
}
