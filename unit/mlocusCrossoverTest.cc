//! \file mlocusCrossoverTest.cc \ingroup unit

#define BOOST_TEST_MODULE mlocusCrossoverTest
#define BOOST_TEST_DYN_LINK 

#include <config.h>
#include <fwdpp/diploid.hh>
#include <boost/test/unit_test.hpp>
#include <unistd.h>
#include <iterator>
#include <functional>
#include <list>

using mut = KTfwd::mutation;
using gtype = KTfwd::gamete;
