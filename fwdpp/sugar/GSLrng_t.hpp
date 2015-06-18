/*! \file GSLrng_t.hpp
  \brief Wrapper for gsl_rng *
  \ingroup sugar
 */
#ifndef __FWDPP_SUGAR_GSLRNG_T_HPP__
#define __FWDPP_SUGAR_GSLRNG_T_HPP__

#include <fwdpp/sugar/gsl/tags.hpp>
#include <fwdpp/sugar/gsl/deleter.hpp>

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
    template<typename U=T>
    sugar::gsl_rng_ptr_t setup(typename std::enable_if<std::is_same<U,GSL_RNG_TAUS2>::value,U>::type)
    {
      return sugar::gsl_rng_ptr_t(gsl_rng_alloc(gsl_rng_taus2));
    }
    template<typename U=T>
    sugar::gsl_rng_ptr_t setup(typename std::enable_if<std::is_same<U,GSL_RNG_MT19937>::value,U>::type)
    {
      return sugar::gsl_rng_ptr_t(gsl_rng_alloc(gsl_rng_mt19937));
    }
  public:
    //! Smart pointer wrapping the gsl_rng *
    sugar::gsl_rng_ptr_t r;

    //! Construct with a seed
    GSLrng_t(const unsigned & seed) : r(setup(T())) {
      gsl_rng_set(r.get(),seed);
    }

    //! Copy constructor
    GSLrng_t( const GSLrng_t & __rng) : r(gsl_rng_clone(__rng.r.get()))
    {
    }

    //! Return underlying pointer
    gsl_rng * get() const {
      return r.get();
    }
  };
}

#endif
