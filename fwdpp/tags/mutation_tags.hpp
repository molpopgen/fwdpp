#ifndef __FWDPP_TAGS_MUTATION_TAGS_HPP__
#define __FWDPP_TAGS_MUTATION_TAGS_HPP__

namespace KTfwd {
  namespace tags {
    //Dispatch tags
    template<bool> struct mmodel_type{};
    using gamete_independent = mmodel_type<false>;
    using gamete_dependent = mmodel_type<true>;
  }
}

#endif
