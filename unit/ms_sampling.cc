/*! 
  \file ms_sampling.cc 
  \ingroup unit 
  \brief Testing KTfwd::ms_sample
*/
#define BOOST_TEST_MODULE ms_sampling
#define BOOST_TEST_DYN_LINK 

#include <iostream>
#include <config.h>
#include <boost/test/unit_test.hpp>
#include <fwdpp/diploid.hh>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/sugar/singlepop.hpp>
#include <fwdpp/sugar/infsites.hpp>
#include <fwdpp/sampling_functions.hpp>

using mutation_t = KTfwd::popgenmut;
using poptype = KTfwd::singlepop<mutation_t>;

void add_mutations(poptype & p)
{
  //auto first_gam = p.gametes.begin();
  for( unsigned i = 0 ; i < 10 ; ++i )
    {
      p.mutations.emplace_back(double(i),0.,0.,1);
      poptype::gamete_t ng(1);
      ng.mutations.push_back(p.mutations.size()-1);

      p.mutations.emplace_back(double(i)+0.5,0.,0.,1);
      poptype::gamete_t ng2(1);
      ng2.mutations.push_back(p.mutations.size()-1);

      p.gametes.emplace_back(std::move(ng));
      p.diploids[i].first = p.gametes.size()-1;
      p.gametes.emplace_back(std::move(ng2));
      p.diploids[i].second = p.gametes.size()-1;
    }
  p.gametes[0].n=0;
  p.mcounts = decltype(p.mcounts)(p.mcounts.size(),1);
  //p.gametes.erase(first_gam);
}

poptype make_poptype()
{
  poptype pop(10);
  add_mutations(pop);
  return pop;
}

//This is making sure that the pop is setup ok.
//Essentially an API check
BOOST_AUTO_TEST_CASE( check_pop )
{
  poptype pop = make_poptype();
  
  BOOST_REQUIRE_EQUAL(pop.mutations.size(),20);
  BOOST_REQUIRE_EQUAL(pop.gametes.size(),21);
}

BOOST_AUTO_TEST_CASE( sample_properties_test )
{
  poptype pop = make_poptype();
  auto sample = KTfwd::fwdpp_internal::ms_sample_separate_single_deme(pop.mutations,pop.gametes,pop.diploids,
								      std::vector<unsigned>({0,1,2,3,4,5,6,7,8,9}),
								      20,true);
  BOOST_REQUIRE_EQUAL(sample.second.empty(),true); //no selected mutations
  BOOST_REQUIRE_EQUAL( std::is_sorted(sample.first.begin(),
				      sample.first.end(),
				      [](const std::pair<double,std::string> & a,
					 const std::pair<double,std::string> & b) {
					return a.first < b.first;
				      } ), true );
  for(unsigned i=0;i<sample.first.size();++i)
    {
      BOOST_REQUIRE_EQUAL(sample.first[i].second.size(),20);
    }
  BOOST_REQUIRE_EQUAL(sample.first.size(),20); //20 neutral mutations
}


BOOST_AUTO_TEST_CASE( odd )
{
  poptype pop = make_poptype();
  auto sample = KTfwd::fwdpp_internal::ms_sample_separate_single_deme(pop.mutations,pop.gametes,pop.diploids,
								      std::vector<unsigned>({0,1,2,3,4,5,6,7,8,9}),
								      19,true);
  BOOST_REQUIRE_EQUAL(sample.second.empty(),true); //no selected mutations
  BOOST_REQUIRE_EQUAL( std::is_sorted(sample.first.begin(),
				      sample.first.end(),
				      [](const std::pair<double,std::string> & a,
					 const std::pair<double,std::string> & b) {
					return a.first < b.first;
				      } ), true );
  for(unsigned i=0;i<sample.first.size();++i)
    {
      BOOST_REQUIRE_EQUAL(sample.first[i].second.size(),19);
    }
  //Aha--need to have a "remove invariant"
  BOOST_REQUIRE_EQUAL(sample.first.size(),19); //20 neutral mutations
}
