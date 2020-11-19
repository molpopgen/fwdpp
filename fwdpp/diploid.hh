/*! \file diploid.hh
  \brief Main header for programming using this library
  \warning Do not try to include individual headers a la carte. Life will get
  confusing for you. Just include this file.

  \code
  #include <fwdpp/diploid.hh>
  \endcode
 */
#ifndef __DIPLOID_HH__
#define __DIPLOID_HH__

#warning("This header is deprecated and will be removed by 0.10.0 or soon thereafter")

// namespace std
#include <vector>
#include <list>
#include <iterator>
#include <ctime>
#include <cmath>

// gsl
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// headers from this project
#include <fwdpp/type_traits.hpp>
#include <fwdpp/debug.hpp>
#include <fwdpp/fitness_models.hpp>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/fwd_functional.hpp>
#include <fwdpp/insertion_policies.hpp>
#include <fwdpp/sampling_functions.hpp>
#include <fwdpp/util.hpp>
#include <fwdpp/sample_diploid.hpp>
#endif

/*! \namespace fwdpp
  \brief The primary namespace defined by this library.
 */

/*! \namespace fwdpp::fwdpp_internal
  \brief Nested namespace for nuts and bolts of certain library functions
*/

/*! \namespace fwdpp::traits
  \brief Nested namespace type traits
*/

/*! \namespace fwdpp::traits::internal
  \brief Nested namespace implementation details of type traits
*/

/*! \namespace fwdpp::tags
  \brief Nested namespace for dispatch tags for template functions.
*/

/*! @defgroup sugar Syntactic sugar layer
  \brief Syntactic sugar for easier development of simulations

  See @ref md_md_sugar for a full description of the features that fwdpp's
  sugar layer provides.
*/
