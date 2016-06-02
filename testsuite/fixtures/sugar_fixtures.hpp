#ifndef FWDPP_TESTSUITE_SUGAR_FIXTURES_HPP
#define FWDPP_TESTSUITE_SUGAR_FIXTURES_HPP

#include <fwdpp/sugar/popgenmut.hpp>
#include <fwdpp/sugar/singlepop.hpp>
#include <fwdpp/sugar/metapop.hpp>

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

#endif
