#ifndef __FWDPP_SUGAR_METAPOP_HPP__
#define __FWDPP_SUGAR_METAPOP_HPP__

#include <fwdpp/sugar/metapop/metapop.hpp>
#include <fwdpp/fwd_functional.hpp>

#ifdef FWDPP_SUGAR_USE_BOOST
#include <boost/container/list.hpp>
#include <boost/container/vector.hpp>
#include <boost/pool/pool_alloc.hpp>
#include <boost/unordered_set.hpp>
#include <boost/functional/hash.hpp>

namespace KTfwd
{
  template<typename mtype> using metapop_mlist_t = boost::container::list<mtype,boost::pool_allocator<mtype> >;
  template<typename mtype> using metapop_gamete_t = gamete_base<mtype,metapop_mlist_t<mtype>>;
  template<typename mtype> using metapop_glist_t = boost::container::list<metapop_gamete_t<mtype>, boost::pool_allocator<metapop_gamete_t<mtype>>>;
  template<typename mtype> using metapop_dipvector_t = boost::container::vector<std::pair<typename metapop_glist_t<mtype>::iterator,
											  typename metapop_glist_t<mtype>::iterator> >;
  template<typename mtype> using metapop = sugar::metapop<mtype,
							  metapop_mlist_t<mtype>,
							  metapop_glist_t<mtype>,
							  metapop_dipvector_t<mtype>,
							  boost::container::vector<metapop_glist_t<mtype>>,
							  boost::container::vector<metapop_dipvector_t<mtype>>,
							  boost::container::vector<mtype>,
							  boost::container::vector<unsigned>,
							  boost::unordered_set<double,boost::hash<double>,KTfwd::equal_eps>>;
}
#else

#include <list>
#include <vector>
#include <unordered_set>
#include <fwdpp/fwd_functional.hpp> 
#include <fwdpp/sugar/metapop/metapop.hpp>
namespace KTfwd
{
  template<typename mtype> using metapop_mlist_t = std::list<mtype>;
  template<typename mtype> using metapop_gamete_t = gamete_base<mtype,metapop_mlist_t<mtype>>;
  template<typename mtype> using metapop_glist_t = std::list<metapop_gamete_t<mtype> >;
  template<typename mtype> using metapop_dipvector_t = std::vector<std::pair<typename metapop_glist_t<mtype>::iterator,
									     typename metapop_glist_t<mtype>::iterator> >;
  template<typename mtype> using metapop = sugar::metapop<mtype,
							  metapop_mlist_t<mtype>,
							  metapop_glist_t<mtype>,
							  metapop_dipvector_t<mtype>,
							  std::vector<metapop_glist_t<mtype>>,
							  std::vector<metapop_dipvector_t<mtype>>,
							  std::vector<mtype>,
							  std::vector<unsigned>,
							  std::unordered_set<double,std::hash<double>,KTfwd::equal_eps>>;
}
#endif

#endif
