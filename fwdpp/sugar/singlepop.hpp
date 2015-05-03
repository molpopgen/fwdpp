#ifndef __FWDPP_SUGAR_SINGLEPOP_HPP__
#define __FWDPP_SUGAR_SINGLEPOP_HPP__

#include <fwdpp/sugar/singlepop/singlepop.hpp>
#include <fwdpp/fwd_functional.hpp>

#ifdef FWDPP_SUGAR_USE_BOOST
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

  template<typename T>
  struct diploid_hack
  {
    T first,second;
    diploid_hack() : first(T()),second(T()) {}
    diploid_hack(  T const & f, T const & s ) : first(f),second(s)
    {
    }
    diploid_hack( const diploid_hack & ) = default;
    diploid_hack( diploid_hack && ) = default;
    diploid_hack & operator=( diploid_hack && ) = default;
    diploid_hack & operator=( diploid_hack const & ) = default;
  };
  template<typename mtype> using DIPT = diploid_hack<typename glist_t<mtype>::iterator>;
  //template<typename mtype> using DIPT = std::pair<typename glist_t<mtype>::iterator,typename glist_t<mtype>::iterator>;
  /*!
    \brief Single locus, single population without serialization.  Cannot be copied, etc.
    See @ref md_md_sugar for rationale, etc.
    \ingroup sugar
  */
  boost::container::vector<DIPT<mtype>> foobar;
  template<typename mtype> using singlepop = sugar::singlepop<mtype,
							      mlist_t<mtype>,
							      glist_t<mtype>,
							      boost::container::vector<DIPT<mtype>>,
							      //boost::container::vector<std::pair<typename glist_t<mtype>::iterator,
							      //typename glist_t<mtype>::iterator> >,
							      boost::container::vector<mtype>,
							      boost::container::vector<unsigned>,
							      boost::unordered_set<double,boost::hash<double>,KTfwd::equal_eps>
							      >;

  /*!
    \brief Single locus, single population with serialization.  Can be copied, etc.
    See @ref md_md_sugar for rationale, etc.
    \ingroup sugar
  */
  template<typename mtype,
	   typename mwriter_t,
	   typename mreader_t> using singlepop_serialized = sugar::singlepop_serialized<mtype,
											mwriter_t,mreader_t,
											mlist_t<mtype>,
											glist_t<mtype>,
											boost::container::vector<std::pair<typename glist_t<mtype>::iterator,
															   typename glist_t<mtype>::iterator> >,
											boost::container::vector<mtype>,
											boost::container::vector<unsigned>,
											boost::unordered_set<double,boost::hash<double>,KTfwd::equal_eps>
											>;
}
#else

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
  /*!
    \brief Single locus, single population without serialization.  Cannot be copied, etc.
    See @ref md_md_sugar for rationale, etc.
    \ingroup sugar
  */
  template<typename mtype> using singlepop = sugar::singlepop<mtype,
							      mlist_t<mtype>,
							      glist_t<mtype>,
							      std::vector< std::pair<typename glist_t<mtype>::iterator,
										     typename glist_t<mtype>::iterator> >,
							      std::vector<mtype>,
							      std::vector<unsigned>,
							      std::unordered_set<double,std::hash<double>,KTfwd::equal_eps>
							      >;

  /*! 
    \brief Single locus, single population with serialization.  Can be copied, etc.
    See @ref md_md_sugar for rationale, etc.
    \ingroup sugar
  */
  template<typename mtype,
	   typename mwriter_t,
	   typename mreader_t> using singlepop_serialized = sugar::singlepop_serialized<mtype,
											mwriter_t,mreader_t,
											mlist_t<mtype>,
											glist_t<mtype>,
											std::vector< std::pair<typename glist_t<mtype>::iterator,
													       typename glist_t<mtype>::iterator> >,
											std::vector<mtype>,
											std::vector<unsigned>,
											std::unordered_set<double,std::hash<double>,KTfwd::equal_eps>
											>;
}
#endif

#endif

