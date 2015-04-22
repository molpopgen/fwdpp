#ifndef __FWDPP_SUGAR_GSL_TAGS_HPP__
#define __FWDPP_SUGAR_GSL_TAGS_HPP__

namespace KTfwd {
  namespace sugar {
    enum class GSL_RNG_TYPE { MT19937, TAUS2 };
    template<GSL_RNG_TYPE> struct GSL_RNG_TYPE_TAG {};
  }
}

#endif
