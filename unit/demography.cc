

/*! \file demography.cc
  \ingroup unit 
  \brief Testing functions in fwdpp/demography.hpp
*/
#define BOOST_TEST_MODULE demography
#define BOOST_TEST_DYN_LINK 

#include <config.h>
#include <algorithm>
#include <boost/test/unit_test.hpp>
#include <fwdpp/diploid.hh>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/sugar/metapop.hpp>
#include <fwdpp/sugar/singlepop.hpp>
#include <fwdpp/sugar/infsites.hpp>
#include <fwdpp/sugar/serialization.hpp>

using mtype = KTfwd::popgenmut;
using metapop_t = KTfwd::metapop<mtype>;

/*
  The first unit tests are of the low-level functions.
  We use a metapop_t for convenience in initializing all the
  required containers, etc.  HOWEVER, these low-level functions
  leave some of the data in a metapop_t inconsistent with the actual data.
  Specifically, metapop_t::Ns is not congruent with the sizes of elements
  in metapop_t::diploids
*/

BOOST_AUTO_TEST_CASE( low_level_copy_deme_test )
{
  {
    metapop_t mpop({1000});
    BOOST_REQUIRE(KTfwd::check_sum(mpop.gametes,2000));
    auto rv = KTfwd::copy_deme(mpop.mutations,mpop.mcounts,mpop.gametes,mpop.diploids,0); //append a copy of deme 0 to the metapop
    BOOST_REQUIRE_EQUAL(rv,0);
    BOOST_REQUIRE(KTfwd::check_sum(mpop.gametes,4000));
    BOOST_REQUIRE_EQUAL(mpop.diploids.size(),2);
    BOOST_REQUIRE(mpop.diploids[0]==mpop.diploids[1]);
    for( const auto & dips : mpop.diploids )
      {
	BOOST_REQUIRE( KTfwd::popdata_sane(dips,mpop.gametes,mpop.mcounts) );
      }
  }
  
  //indexes out of range
  {
    metapop_t mpop({1000});
    metapop_t mpop_copy(mpop); //make a copy to check that return on error leaves mpop intact
    auto rv = KTfwd::copy_deme(mpop.mutations,mpop.mcounts,mpop.gametes,mpop.diploids,1); //Now, 1 is out of range, so we'll get an error
    BOOST_REQUIRE(rv != 0);
    BOOST_REQUIRE(mpop == mpop_copy);
  }
}

BOOST_AUTO_TEST_CASE( low_level_merge_deme_test )
{
  //Two-deme case
  {
    metapop_t mpop{1000,500};
    auto rv = KTfwd::merge_demes(mpop.diploids,0,1);
    BOOST_REQUIRE_EQUAL(rv,0);
    BOOST_REQUIRE_EQUAL(mpop.diploids.size(),1);
    BOOST_REQUIRE_EQUAL(mpop.diploids[0].size(),1500);
    BOOST_REQUIRE(KTfwd::check_sum(mpop.gametes,3000));
  }

  //Three-deme case
  {
    metapop_t mpop{1000,500,250};
    auto rv = KTfwd::merge_demes(mpop.diploids,1,2);
    BOOST_REQUIRE_EQUAL(rv,0);
    BOOST_REQUIRE_EQUAL(mpop.diploids.size(),2);
    BOOST_REQUIRE_EQUAL(mpop.diploids[1].size(),750);
    BOOST_REQUIRE(KTfwd::check_sum(mpop.gametes,3500));
  }

  //Will give same result as previous test because 2 and 1 will get swapped
  {
    metapop_t mpop{1000,500,250};
    auto rv = KTfwd::merge_demes(mpop.diploids,2,1);
    BOOST_REQUIRE_EQUAL(rv,0);
    BOOST_REQUIRE_EQUAL(mpop.diploids.size(),2);
    BOOST_REQUIRE_EQUAL(mpop.diploids[1].size(),750);
    BOOST_REQUIRE(KTfwd::check_sum(mpop.gametes,3500));
  }

  //Will not do anything b/c i == j
  {
    metapop_t mpop{1000,500,250};
    metapop_t mpop_copy(mpop);
    auto rv = KTfwd::merge_demes(mpop.diploids,1,1);
    BOOST_REQUIRE_EQUAL(rv,1);
    BOOST_REQUIRE(mpop==mpop_copy);
  } 

  //Two-deme case with index out of range
  {
    metapop_t mpop{1000,500};
    metapop_t mpop_copy(mpop);
    auto rv = KTfwd::merge_demes(mpop.diploids,0,2); //2 is out of range
    BOOST_REQUIRE_EQUAL(rv,-1);
    BOOST_REQUIRE(mpop==mpop_copy);
  }
}

BOOST_AUTO_TEST_CASE( low_level_remove_deme_test )
{
  {
    metapop_t mpop{1000,2000,3000};
    //remove the second deme:
    auto rv = KTfwd::remove_deme(mpop.mutations,mpop.mcounts,mpop.gametes,mpop.diploids,1);
    BOOST_REQUIRE(rv==0);
    BOOST_REQUIRE(mpop.diploids.size()==2);
    BOOST_REQUIRE(mpop.diploids[0].size()==1000);
    BOOST_REQUIRE(mpop.diploids[1].size()==3000);
    BOOST_REQUIRE(KTfwd::check_sum(mpop.gametes,8000));
  }
  
  {
    metapop_t mpop{1000,2000,3000};
    metapop_t mpop_copy(mpop);
    //deme index out of range:
    auto rv = KTfwd::remove_deme(mpop.mutations,mpop.mcounts,mpop.gametes,mpop.diploids,3);
    BOOST_REQUIRE(rv==-1);
    BOOST_REQUIRE(KTfwd::check_sum(mpop.gametes,12000));
    BOOST_REQUIRE(mpop==mpop_copy);
  }
}

BOOST_AUTO_TEST_CASE( low_level_swap_demes_test )
{
  {
    metapop_t mpop{1000,2000,3000,4000,5000};
    auto rv = KTfwd::swap_demes(mpop.diploids,0,4);
    BOOST_REQUIRE_EQUAL(rv,0);
    BOOST_REQUIRE(KTfwd::check_sum(mpop.gametes,30000));
    BOOST_REQUIRE_EQUAL(mpop.diploids[0].size(),5000);
    BOOST_REQUIRE_EQUAL(mpop.diploids[4].size(),1000);
  }

  {
    metapop_t mpop{1000,2000,3000,4000,5000};
    metapop_t mpop_copy(mpop);
    auto rv = KTfwd::swap_demes(mpop.diploids,1,1); //i==j, so nothing will happen
    BOOST_REQUIRE_EQUAL(rv,1);
    BOOST_REQUIRE(KTfwd::check_sum(mpop.gametes,30000));
    BOOST_REQUIRE(mpop==mpop_copy);
  }
  
  {
    metapop_t mpop{1000,2000,3000,4000,5000};
    metapop_t mpop_copy(mpop);
    auto rv = KTfwd::swap_demes(mpop.diploids,0,5); //deme index j out of range
    BOOST_REQUIRE_EQUAL(rv,-1);
    BOOST_REQUIRE(KTfwd::check_sum(mpop.gametes,30000));
    BOOST_REQUIRE(mpop==mpop_copy);
  }
}

BOOST_AUTO_TEST_CASE( low_level_split_demes_test )
{
  {
    metapop_t mpop{1000};
    KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);
    auto rv = KTfwd::split_deme(rng.get(),mpop.mutations,mpop.mcounts,mpop.gametes,mpop.diploids,0,250,false);
    BOOST_REQUIRE_EQUAL(rv,0);
    BOOST_REQUIRE_EQUAL(mpop.diploids.size(),2);
    BOOST_REQUIRE_EQUAL(mpop.diploids[0].size(),750);
    BOOST_REQUIRE_EQUAL(mpop.diploids[1].size(),250);
    BOOST_REQUIRE(KTfwd::check_sum(mpop.gametes,2000));
  }

  //with replacement
  {
    metapop_t mpop{1000};
    KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);
    auto rv = KTfwd::split_deme(rng.get(),mpop.mutations,mpop.mcounts,mpop.gametes,mpop.diploids,0,990,true);
    BOOST_REQUIRE_EQUAL(rv,0);
    BOOST_REQUIRE_EQUAL(mpop.diploids.size(),2);
    BOOST_REQUIRE_EQUAL(mpop.diploids[0].size(),10);
    BOOST_REQUIRE_EQUAL(mpop.diploids[1].size(),990);
    BOOST_REQUIRE(KTfwd::check_sum(mpop.gametes,2000));
  }

  //index out of range
  {
    metapop_t mpop{1000};
    metapop_t mpop_copy(mpop);
    KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);
    auto rv = KTfwd::split_deme(rng.get(),mpop.mutations,mpop.mcounts,mpop.gametes,mpop.diploids,1,250,false);
    BOOST_REQUIRE_EQUAL(rv,-1);
    BOOST_REQUIRE(mpop==mpop_copy);
  }
}

BOOST_AUTO_TEST_CASE( low_level_admix_demes_test )
{
  {
    metapop_t mpop{1000,1000};
    KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);
    auto rv = KTfwd::admix_demes(rng.get(),mpop.mutations,mpop.mcounts,mpop.gametes,mpop.diploids,0,1,0.25,1000,false);
    BOOST_REQUIRE_EQUAL(rv,0);
    BOOST_REQUIRE(KTfwd::check_sum(mpop.gametes,6000));
  }

  //with resampling
  {
    metapop_t mpop{1000,1000};
    KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);
    auto rv = KTfwd::admix_demes(rng.get(),mpop.mutations,mpop.mcounts,mpop.gametes,mpop.diploids,0,1,0.25,1000,true);
    BOOST_REQUIRE_EQUAL(rv,0);
    BOOST_REQUIRE(KTfwd::check_sum(mpop.gametes,6000));
  }

  //a pop index out of range
  {
    metapop_t mpop{1000,1000};
    KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);
    auto rv = KTfwd::admix_demes(rng.get(),mpop.mutations,mpop.mcounts,mpop.gametes,mpop.diploids,0,2,0.25,1000,false);
    BOOST_REQUIRE_EQUAL(rv,-1);
  }

  //Fraction of pop i ancestry < 0.
  {
    metapop_t mpop{1000,1000};
    KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);
    auto rv = KTfwd::admix_demes(rng.get(),mpop.mutations,mpop.mcounts,mpop.gametes,mpop.diploids,0,1,-0.25,1000,false);
    BOOST_REQUIRE_EQUAL(rv,1);
  }

  //Fraction of pop i ancestry >= 1.
  {
    metapop_t mpop{1000,1000};
    KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);
    auto rv = KTfwd::admix_demes(rng.get(),mpop.mutations,mpop.mcounts,mpop.gametes,mpop.diploids,0,1,1.,1000,false);
    BOOST_REQUIRE_EQUAL(rv,1);
  }

  /*
    When sampling w/o replacement, the new population size cannot be so large
    that the number of individuals required from deme i or j is >= than those demes'
    respective sizes.  This test samples so many individuals that this condition cannot be met,
    resulting in a return value of 1
  */

  {
    metapop_t mpop{1000,1000};
    KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);
    //here, we want a new deme that is 50% ancestry from demes 0 and 1.  But, sampling w/o
    //replacement for 2000 individuals will require all individuals from demes i and j.
    auto rv = KTfwd::admix_demes(rng.get(),mpop.mutations,mpop.mcounts,mpop.gametes,mpop.diploids,0,1,0.5,2000,false);
    BOOST_REQUIRE_EQUAL(rv,1);
  }

  //sampling w/replacement has no such limitations
  {
    metapop_t mpop{1000,1000};
    KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);
    //here, we want a new deme that is 50% ancestry from demes 0 and 1.  But, sampling w/o
    //replacement for 2000 individuals will require all individuals from demes i and j.
    auto rv = KTfwd::admix_demes(rng.get(),mpop.mutations,mpop.mcounts,mpop.gametes,mpop.diploids,0,1,0.5,2000,true);
    BOOST_REQUIRE_EQUAL(rv,0);
  }
}

#include <fwdpp/sugar/demography.hpp>

/*
  Now, we test the higher-level functions in the sugar sub-library.

  These functions update other data inside of KTfwd::metapop.

  These functions below are implemented via calls to the above functions.

  Thus, the tests below are API tests plus tests that KTfwd::sugar::metapop::Ns
  is getting properly updated..
*/
BOOST_AUTO_TEST_CASE( copy_pop_test )
{
  {
    metapop_t mpop{1000};
    auto rv = KTfwd::copy_pop(mpop,0);
    BOOST_REQUIRE_EQUAL(rv,0);
    BOOST_REQUIRE_EQUAL(mpop.Ns.size(),2);
    BOOST_REQUIRE_EQUAL(mpop.Ns.size(),mpop.diploids.size());
    for(std::size_t i = 0 ; i < mpop.Ns.size() ; ++i)
      {
	BOOST_REQUIRE_EQUAL(mpop.Ns[i],mpop.diploids[i].size());
      }
  }
}

BOOST_AUTO_TEST_CASE( merge_pops_test )
{
  {
    metapop_t mpop({1000,500});
    auto rv = KTfwd::merge_pops(mpop,0,1);
    BOOST_REQUIRE_EQUAL(rv,0);
    BOOST_REQUIRE(KTfwd::check_sum(mpop.gametes,3000));
    BOOST_REQUIRE_EQUAL(mpop.Ns.size(),1);
    BOOST_REQUIRE_EQUAL(mpop.Ns.size(),mpop.diploids.size());
    for(std::size_t i = 0 ; i < mpop.Ns.size() ; ++i)
      {
	BOOST_REQUIRE_EQUAL(mpop.Ns[i],mpop.diploids[i].size());
      }
  }
}

BOOST_AUTO_TEST_CASE( remove_pop_test )
{
  {
    metapop_t mpop({1000,500});
    auto rv = KTfwd::remove_pop(mpop,1);
    BOOST_REQUIRE_EQUAL(rv,0);
    BOOST_REQUIRE(KTfwd::check_sum(mpop.gametes,2000));
    BOOST_REQUIRE_EQUAL(mpop.Ns.size(),1);
    BOOST_REQUIRE_EQUAL(mpop.Ns.size(),mpop.diploids.size());
    for(std::size_t i = 0 ; i < mpop.Ns.size() ; ++i)
      {
	BOOST_REQUIRE_EQUAL(mpop.Ns[i],mpop.diploids[i].size());
      }
  }
}

BOOST_AUTO_TEST_CASE( split_pop_test )
{
  {
    metapop_t mpop{1000};
    KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);
    auto rv = KTfwd::split_pop(rng.get(),mpop,0,500,false);
    BOOST_REQUIRE_EQUAL(rv,0);
    BOOST_REQUIRE(KTfwd::check_sum(mpop.gametes,2000));
    BOOST_REQUIRE_EQUAL(mpop.Ns.size(),2);
    BOOST_REQUIRE_EQUAL(mpop.Ns.size(),mpop.diploids.size());
    for(std::size_t i = 0 ; i < mpop.Ns.size() ; ++i)
      {
	BOOST_REQUIRE_EQUAL(mpop.Ns[i],mpop.diploids[i].size());
      }
  }
  {
    metapop_t mpop{1000};
    KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);
    auto rv = KTfwd::split_pop(rng.get(),mpop,0,500,true);
    BOOST_REQUIRE_EQUAL(rv,0);
    BOOST_REQUIRE(KTfwd::check_sum(mpop.gametes,2000));
    BOOST_REQUIRE_EQUAL(mpop.Ns.size(),2);
    BOOST_REQUIRE_EQUAL(mpop.Ns.size(),mpop.diploids.size());
    for(std::size_t i = 0 ; i < mpop.Ns.size() ; ++i)
      {
	BOOST_REQUIRE_EQUAL(mpop.Ns[i],mpop.diploids[i].size());
      }
  }
}

BOOST_AUTO_TEST_CASE( admix_pops_test )
{
  {
    metapop_t mpop{1000,1000};
    KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);
    auto rv = KTfwd::admix_pops(rng.get(),mpop,0,1,0.5,1000);
    BOOST_REQUIRE_EQUAL(rv,0);
    BOOST_REQUIRE(KTfwd::check_sum(mpop.gametes,6000));
    BOOST_REQUIRE_EQUAL(mpop.Ns.size(),3);
    BOOST_REQUIRE_EQUAL(mpop.Ns.size(),mpop.diploids.size());
    for(std::size_t i = 0 ; i < mpop.Ns.size() ; ++i)
      {
	BOOST_REQUIRE_EQUAL(mpop.Ns[i],mpop.diploids[i].size());
      }
  }
}
