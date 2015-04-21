#ifndef __FWDPP_SUGAR_SINGLEPOP_BOOST_HPP__
#define __FWDPP_SUGAR_SINGLEPOP_BOOST_HPP__

#include <fwdpp/sugar/singlepop/singlepop.hpp>
#include <fwdpp/fwd_functional.hpp> 
#include <boost/container/list.hpp>
#include <boost/container/vector.hpp>
#include <boost/pool/pool_alloc.hpp>
#include <boost/unordered_set.hpp>
#include <boost/functional/hash.hpp>

namespace KTfwd
{
  template<typename mtype> using mlist_t = boost::container::list<mtype,boost::pool_allocator<mtype> >;
  template<typename mtype> using gamete_t = gamete_base<mtype,mlist_t<mtype>>;
  template<typename mtype> using glist_t = boost::container::list<gamete_t<mtype>, boost::pool_allocator<gamete_t<mtype>>>;
  template<typename mtype> using singlepop = sugar::singlepop<mtype,
							      mlist_t<mtype>,
							      glist_t<mtype>,
							      boost::container::vector<mtype>,
							      boost::container::vector<unsigned>,
							      boost::unordered_set<double,boost::hash<double>,KTfwd::equal_eps>
							      >;

  template<typename mtype,
	   typename mwriter_t,
	   typename mreader_t> using singlepop_serialized = sugar::singlepop_serialized<mtype,
											mwriter_t,mreader_t,
											mlist_t<mtype>,
											glist_t<mtype>,
											boost::container::vector<mtype>,
											boost::container::vector<unsigned>,
											boost::unordered_set<double,boost::hash<double>,KTfwd::equal_eps>
											>;
}
#endif

