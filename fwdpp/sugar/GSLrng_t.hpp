/*! \file GSLrng_t.hpp
  \brief Wrapper for gsl_rng *
  \ingroup sugar
 */
#ifndef __FWDPP_SUGAR_GSLRNG_T_HPP__
#define __FWDPP_SUGAR_GSLRNG_T_HPP__

#include <fwdpp/sugar/gsl/tags.hpp>
#include <fwdpp/sugar/gsl/deleter.hpp>
#include <cassert>

namespace KTfwd {

  //! Distpatch tag to signal GSLrng_t to instantiate in terms of gsl_rng_mt19937
  using GSL_RNG_MT19937 = sugar::GSL_RNG_TYPE_TAG<sugar::GSL_RNG_TYPE::MT19937>;
    //! Distpatch tag to signal GSLrng_t to instantiate in terms of gsl_rng_taus2
  using GSL_RNG_TAUS2 = sugar::GSL_RNG_TYPE_TAG<sugar::GSL_RNG_TYPE::TAUS2>;

  /*!
    \brief A wrapper around gsl_rng * objects.
    
    The template instantiation type must be a model of 
    KTfwd::sugar::GSL_RNG_TYPE_TAG, which specifies the
    gsl_rng type.

    This type holds an object of type KTfwd::sugar::gsl_rng_ptr_t,
    which is a smart pointer that manages freeing the gsl_rng * upon
    destruction.
   */
  template<typename T>
  class GSLrng_t {
  private:
    sugar::gsl_rng_ptr_t setup( GSL_RNG_MT19937 )
    {
      return sugar::gsl_rng_ptr_t(gsl_rng_alloc(gsl_rng_mt19937));
    }

    sugar::gsl_rng_ptr_t setup( GSL_RNG_TAUS2 )
    {
      return sugar::gsl_rng_ptr_t(gsl_rng_alloc(gsl_rng_taus2));
    }
  public:
    //! Smart pointer wrapping the gsl_rng *
    sugar::gsl_rng_ptr_t r;

    //! Typedef for RNG type, if needed
    using rngtype = T;
    
    //! Construct with a seed
    GSLrng_t(const unsigned long seed) : r(setup(T())) {
      gsl_rng_set(r.get(),seed);
    }

    //! Copy constructor
    GSLrng_t( const GSLrng_t & rng) : r(setup(T())) {
      int rv = gsl_rng_memcpy(r.get(),rng.get());
      assert(rv==GSL_SUCCESS);
    }
      
    //! Copy constructor
    GSLrng_t( GSLrng_t & rng) : r(setup(T())) {
      int rv = gsl_rng_memcpy(r.get(),rng.get());
      assert(rv==GSL_SUCCESS);
    }

    GSLrng_t & operator=(GSLrng_t &) = delete;
    
    //! Return underlying pointer
    gsl_rng * get() const {
      return r.get();
    }
  };
}

#endif
