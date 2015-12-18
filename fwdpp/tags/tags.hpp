#ifndef __FWDPP_TAGS_TAGS_HPP__
#define __FWDPP_TAGS_TAGS_HPP__

namespace KTfwd {
  namespace tags {
    //! Dispatch tags for mutations/gametes
    template<bool> struct extinct_t{};
    //! Signals constructor for an extinct type
    using extinct = extinct_t<true>;
    /*! 
      Signals constructor for an extant type. 
      \note This is provided for completeness, 
      but is not expected to be used
    */
    using not_extinct = extinct_t<false>;
  }
}

#endif
