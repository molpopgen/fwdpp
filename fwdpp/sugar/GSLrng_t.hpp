#ifndef __FWDPP_SUGAR_GSLRNG_T_HPP__
#define __FWDPP_SUGAR_GSLRNG_T_HPP__

#include <fwdpp/sugar/gsl/tags.hpp>
#include <fwdpp/sugar/gsl/deleter.hpp>

namespace KTfwd {

  using GSL_RNG_MT19937 = sugar::GSL_RNG_TYPE_TAG<sugar::GSL_RNG_TYPE::MT19937>;
  using GSL_RNG_TAUS2 = sugar::GSL_RNG_TYPE_TAG<sugar::GSL_RNG_TYPE::TAUS2>;
  
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
    sugar::gsl_rng_ptr_t r;

    GSLrng_t(const unsigned & seed) : r(setup(T())) {
      gsl_rng_set(r.get(),seed);
    }
    /*
      boost python requires that things be copy-constructable
      However, the gsl_rng_ptr_t is not, so we gotta write
      a manual copy-constructor
    */
    GSLrng_t( const GSLrng_t & __rng) : r(gsl_rng_clone(__rng.r.get()))
    {
    }
    //Allow explicit conversiont
    operator gsl_rng *() const { return r.get(); }
  };
}

#endif
