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

// namespace std
#include <vector>
#include <list>
#include <cassert>
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
#include <fwdpp/mutation.hpp>
#include <fwdpp/util.hpp>
#include <fwdpp/IO.hpp>
#include <fwdpp/recombination.hpp>
#include <fwdpp/interlocus_recombination.hpp>
#include <fwdpp/sample_diploid.hpp>
#include <fwdpp/demography.hpp>
#endif

/*! \namespace KTfwd
  \brief The primary namespace defined by this library.
 */

/*! \namespace KTfwd::fwdpp_internal
  \brief Nested namespace for nuts and bolts of certain library functions
*/

/*! \namespace KTfwd::traits
  \brief Nested namespace type traits
*/

/*! \namespace KTfwd::traits::internal
  \brief Nested namespace implementation details of type traits
*/

/*! \namespace KTfwd::tags
  \brief Nested namespace for dispatch tags for template functions.
*/

/*! @defgroup sugar Syntactic sugar layer
  \brief Syntactic sugar for easier development of simulations

  See @ref md_md_sugar for a full description of the features that fwdpp's
  sugar layer provides.
 */

/*! @defgroup mlocus Multi-locus/region simulations
 * \brief Functions related to modeling multi-locus/region simulations
 */

/*! \namespace KTfwd::sugar
  \brief Nested namespace for sugar layer.

  This namespace provides the implementation details for @ref sugar.

  See @ref md_md_sugar for a full description of the features that fwdpp's
  sugar layer provides.
 */
