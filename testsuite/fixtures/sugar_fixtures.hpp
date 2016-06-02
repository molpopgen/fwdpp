#ifndef FWDPP_TESTSUITE_SUGAR_FIXTURES_HPP
#define FWDPP_TESTSUITE_SUGAR_FIXTURES_HPP

#include <fwdpp/sugar/popgenmut.hpp>
#include <fwdpp/sugar/singlepop.hpp>
#include <fwdpp/sugar/metapop.hpp>
#include <fwdpp/sugar/multiloc.hpp>

struct singlepop_popgenmut_fixture
{
  KTfwd::singlepop<KTfwd::popgenmut> pop;
  singlepop_popgenmut_fixture() : pop(KTfwd::singlepop<KTfwd::popgenmut>(1000))
  {
  }
};

struct metapop_popgenmut_fixture
{
  KTfwd::metapop<KTfwd::popgenmut> pop;
  metapop_popgenmut_fixture() : pop(KTfwd::metapop<KTfwd::popgenmut>({1000,1000}))
  {
  }
};

struct multiloc_popgenmut_fixture
{
  KTfwd::multiloc<KTfwd::popgenmut> pop;
  multiloc_popgenmut_fixture() : pop(KTfwd::multiloc<KTfwd::popgenmut>(1000,5))
				/*! N=1000, 5 loci */
  {
  }
};

#endif
