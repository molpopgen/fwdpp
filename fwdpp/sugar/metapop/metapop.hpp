#ifndef __FWDPP_SUGAR_METAPOP_HPP__
#define __FWDPP_SUGAR_METAPOP_HPP__

#include <type_traits>
#include <vector>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/sugar/serialization.hpp>

namespace KTfwd {
  namespace sugar {
    /*!
      Abstraction of what is needed to simulate a single population
      using an individual-based sampler from fwdpp
      
      All that is missing is the mutation_type and the container types.
      
      Does not allow copy construction/assignment!
    */
    template<typename mutation_type,
	     typename mlist,
	     typename glist,
	     typename dipvector,
	     typename vglist,
	     typename vdipvector,
	     typename mvector,
	     typename ftvector,
	     typename lookup_table_type>
    class metapop
    {
    public:
      std::vector<unsigned> Ns;

      //Typedefs for various container
      using mutation_t = mutation_type;
      using gamete_t = typename glist::value_type;
      using dipvector_t = dipvector;
      using diploid_t = typename dipvector_t::value_type;
      using mlist_t = mlist;
      using glist_t = glist;
      using vglist_t = vglist;
      using lookup_table_t = lookup_table_type;
  }
}

#endif
