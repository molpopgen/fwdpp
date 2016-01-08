#ifndef __FWDPP_SUGAR_SINGLEPOP_HPP__
#define __FWDPP_SUGAR_SINGLEPOP_HPP__

#include <utility>
#include <vector>
#include <unordered_set>
#include <fwdpp/fwd_functional.hpp>
#include <fwdpp/sugar/singlepop/singlepop.hpp>

namespace KTfwd
{
  template<typename mtype> using singlepop_mvec_t = std::vector<mtype,std::allocator<mtype> >;
  using singlepop_gamete_t = gamete;
  template<typename mtype> using singlepop_gvec_t = std::vector<singlepop_gamete_t, std::allocator<singlepop_gamete_t>>;
  /*!
    \brief Single locus, single population without serialization.  Cannot be copied, etc.
    See @ref md_md_sugar for rationale, etc.
    \ingroup sugar
  */
  template<typename mtype,
	   typename diploid_t = std::pair<std::size_t,std::size_t> >
    using singlepop = sugar::singlepop<mtype,
				       singlepop_mvec_t<mtype>,
				       singlepop_gvec_t<mtype>,
				       std::vector<diploid_t>,
				       std::vector<mtype>,
				       std::vector<uint_t>,
				       std::unordered_set<double,std::hash<double>,KTfwd::equal_eps>
				       >;
}
#endif

