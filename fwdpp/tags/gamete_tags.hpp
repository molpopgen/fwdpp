#ifndef __FWDPP_TAGS_GAMETE_TAGS_HPP__
#define __FWDPP_TAGS_GAMETE_TAGS_HPP__

namespace KTfwd {
  namespace tags {
    //! Built-in types of gametes for which dispatch tags exist
    enum class gamete_type { standard };
    //! Dispatch tags for gamete types
    template<gamete_type> struct gamete_type_tag {};
    //! The "regular" gamete for a diploid simulation
    using standard_gamete = gamete_type_tag<gamete_type::standard>;
  }
}

#endif
