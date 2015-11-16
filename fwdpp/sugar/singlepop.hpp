#ifndef __FWDPP_SUGAR_SINGLEPOP_HPP__
#define __FWDPP_SUGAR_SINGLEPOP_HPP__

#include <utility>
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
  template<typename mtype> using singlepop_mlist_t = boost::container::list<mtype,boost::fast_pool_allocator<mtype> >;
  template<typename mtype> using singlepop_gamete_t = gamete_base<mtype,singlepop_mlist_t<mtype>>;
  template<typename mtype> using singlepop_glist_t = boost::container::list<singlepop_gamete_t<mtype>,
									    boost::fast_pool_allocator<singlepop_gamete_t<mtype>>>;

  /*!
    \brief Single locus, single population without serialization.  Cannot be copied, etc.
    See @ref md_md_sugar for rationale, etc.
    \ingroup sugar
  */
  template<typename mtype,
	   typename diploid_t = std::pair<typename singlepop_glist_t<mtype>::iterator,
					  typename singlepop_glist_t<mtype>::iterator> >
  using singlepop = sugar::singlepop<mtype,
				     singlepop_mlist_t<mtype>,
				     singlepop_glist_t<mtype>,
				     boost::container::vector< diploid_t >,
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
	   typename mreader_t,
	   typename diploid_t = std::pair<typename singlepop_glist_t<mtype>::iterator,
					  typename singlepop_glist_t<mtype>::iterator>,
	   typename diploid_writer_t = diploidIOplaceholder,
	   typename diploid_reader_t = diploidIOplaceholder>
  using singlepop_serialized = sugar::singlepop_serialized<mtype,
							   mwriter_t,mreader_t,
							   singlepop_mlist_t<mtype>,
							   singlepop_glist_t<mtype>,
							   boost::container::vector< diploid_t >,
							   boost::container::vector<mtype>,
							   boost::container::vector<uint_t>,
							   boost::unordered_set<double,boost::hash<double>,KTfwd::equal_eps>,
							   diploid_writer_t,
							   diploid_reader_t
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
  template<typename mtype> using singlepop_mlist_t = std::list<mtype,std::allocator<mtype> >;
  template<typename mtype> using singlepop_gamete_t = gamete_base<mtype,singlepop_mlist_t<mtype>>;
  template<typename mtype> using singlepop_glist_t = std::list<singlepop_gamete_t<mtype>, std::allocator<singlepop_gamete_t<mtype>>>;
  /*!
    \brief Single locus, single population without serialization.  Cannot be copied, etc.
    See @ref md_md_sugar for rationale, etc.
    \ingroup sugar
  */
  template<typename mtype,
	   typename diploid_t = std::pair<typename singlepop_glist_t<mtype>::iterator,
					  typename singlepop_glist_t<mtype>::iterator> >
    using singlepop = sugar::singlepop<mtype,
				       singlepop_mlist_t<mtype>,
				       singlepop_glist_t<mtype>,
				       std::vector< diploid_t >,
				       std::vector<mtype>,
				       std::vector<uint_t>,
				       std::unordered_set<double,std::hash<double>,KTfwd::equal_eps>
				       >;

  /*! 
    \brief Single locus, single population with serialization.  Can be copied, etc.
    See @ref md_md_sugar for rationale, etc.
    \ingroup sugar
  */
  template<typename mtype,
	   typename mwriter_t,
	   typename mreader_t,
	   typename diploid_t = std::pair<typename singlepop_glist_t<mtype>::iterator,
					  typename singlepop_glist_t<mtype>::iterator>,
	   typename diploid_writer_t = diploidIOplaceholder,
	   typename diploid_reader_t = diploidIOplaceholder>
  using singlepop_serialized = sugar::singlepop_serialized<mtype,
							   mwriter_t,mreader_t,
							   singlepop_mlist_t<mtype>,
							   singlepop_glist_t<mtype>,
							   std::vector< diploid_t >,
							   std::vector<mtype>,
							   std::vector<uint_t>,
							   std::unordered_set<double,std::hash<double>,KTfwd::equal_eps>,
							   diploid_writer_t,
							   diploid_reader_t
							   >;
}
#endif

#endif

