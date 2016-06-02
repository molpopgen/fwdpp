#ifndef FWDPP_TESTSUITE_SUGAR_FIXTURES_HPP
#define FWDPP_TESTSUITE_SUGAR_FIXTURES_HPP

#include <fwdpp/sugar/popgenmut.hpp>
#include <fwdpp/sugar/singlepop.hpp>
#include <fwdpp/sugar/metapop.hpp>
#include <fwdpp/sugar/multiloc.hpp>

#include <testsuite/util/custom_dip.hpp>

struct singlepop_popgenmut_fixture
{
  using poptype= KTfwd::singlepop<KTfwd::popgenmut>;
  poptype pop;
  singlepop_popgenmut_fixture() : pop(poptype(1000))
  {
  }
};

struct singlepop_popgenmut_custom_fixture
{
  using poptype = KTfwd::singlepop<KTfwd::popgenmut,custom_diploid_testing_t>;
  poptype pop;
  singlepop_popgenmut_custom_fixture() : pop(poptype(1000))
  {
  }
};

struct metapop_popgenmut_fixture
{
  using poptype = KTfwd::metapop<KTfwd::popgenmut>;
  poptype pop;
  metapop_popgenmut_fixture() : pop(poptype{1000,1000})
  {
  }
};

struct multiloc_popgenmut_fixture
{
  using poptype = KTfwd::multiloc<KTfwd::popgenmut>;
  poptype pop;
  multiloc_popgenmut_fixture() : pop(poptype(1000,5))
				/*! N=1000, 5 loci */
  {
  }
};

#endif
