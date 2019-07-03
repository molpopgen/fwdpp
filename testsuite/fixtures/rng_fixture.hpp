#ifndef FWDPP_TESTSUITE_RNG_FIXTURE_HPP__
#define FWDPP_TESTSUITE_RNG_FIXTURE_HPP__

#include <fwdpp/GSLrng_t.hpp>

struct rng_fixture
{
	fwdpp::GSLrng_mt rng;
	rng_fixture() : rng(42)
	{
	}
};

#endif
