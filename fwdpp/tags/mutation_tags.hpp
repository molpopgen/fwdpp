#ifndef __FWDPP_TAGS_MUTATION_TAGS_HPP__
#define __FWDPP_TAGS_MUTATION_TAGS_HPP__

namespace KTfwd {
  namespace tags {
    //! Dispatch tags for mutation models
    template<bool> struct mmodel_type{};
    //! Gamete-independent mutation model tag (will typically not be used)
    using gamete_independent = mmodel_type<false>;
    //! Gamete-dependent mutation model tag \example HOCind.cc
    using gamete_dependent = mmodel_type<true>;
  }
}

#endif
