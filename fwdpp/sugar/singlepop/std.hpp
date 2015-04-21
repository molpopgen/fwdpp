#ifndef __FWDPP_SUGAR_SINGLEPOP_STD_HPP__
#define __FWDPP_SUGAR_SINGLEPOP_STD_HPP__

#include <list>
#include <vector>
#include <unordered_set>
#include <fwdpp/fwd_functional.hpp> 
#include <fwdpp/sugar/singlepop/singlepop.hpp>

namespace KTfwd
{
  template<typename mtype> using mlist_t = std::list<mtype,std::allocator<mtype> >;
  template<typename mtype> using gamete_t = gamete_base<mtype,mlist_t<mtype>>;
  template<typename mtype> using glist_t = std::list<gamete_t<mtype>, std::allocator<gamete_t<mtype>>>;
  template<typename mtype> using singlepop = sugar::singlepop<mtype,
							      mlist_t<mtype>,
							      glist_t<mtype>,
							      std::vector<mtype>,
							      std::vector<unsigned>,
							      std::unordered_set<double,std::hash<double>,KTfwd::equal_eps>
							      >;
}
#endif

