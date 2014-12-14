/*! \file diploid.hh
  \brief Main header for programming using this library
  \warning Do not try to include individual headers a la carte. Life will get confusing for you. Just include this file.

  \code
  #include <fwdpp/diploid.hh>
  \endcode
 */
#ifndef __DIPLOID_HH__
#define __DIPLOID_HH__

//namespace std
#include <vector>
#include <list>
#include <cassert>
#include <iterator>
#include <ctime>
#include <cmath>

//gsl
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

//headers from this project
#include <fwdpp/debug.hpp>
#include <fwdpp/migration.hpp>
#include <fwdpp/diploid_functions.hpp>
#include <fwdpp/fitness_models.hpp>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/fwd_functional.hpp>
#include <fwdpp/insertion_policies.hpp>
#include <fwdpp/sampling_functions.hpp>
#include <fwdpp/mutation.hpp>
#include <fwdpp/util.hpp>
#include <fwdpp/initms.hpp>
#include <fwdpp/IO.hpp>
#include <fwdpp/recombination_methods_ind.hpp>
#endif

/*! \namespace KTfwd
  \brief The primary namespace defined by this library.
 */

/*! \namespace KTfwd::fwdpp_internal
  \brief Nested namespace for nuts and bolts of certain library functions
*/
