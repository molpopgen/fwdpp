#ifndef __FWDPP_SUGAR_MULTILOC_HPP__
#define __FWDPP_SUGAR_MULTILOC_HPP__

#include <fwdpp/sugar/multiloc/multiloc.hpp>
#include <fwdpp/fwd_functional.hpp>

#ifdef FWDPP_SUGAR_USE_BOOST
#include <boost/container/list.hpp>
#include <boost/container/vector.hpp>
#include <boost/pool/pool_alloc.hpp>
#include <boost/unordered_set.hpp>
#include <boost/functional/hash.hpp>

namespace KTfwd
{
  template<typename mtype> using multiloc_mlist_t = boost::container::list<mtype,boost::fast_pool_allocator<mtype> >;
  template<typename mtype> using multiloc_gamete_t = gamete_base<mtype,multiloc_mlist_t<mtype>>;
  //IDEA:
  template<typename mtype> using multiloc_glist_t = boost::container::list<multiloc_gamete_t<mtype>, boost::fast_pool_allocator<multiloc_gamete_t<mtype>>>;

  /*!
    \brief Single population, multilocus simulation without serialization.  Cannot be copied, etc.
    See @ref md_md_sugar for rationale, etc.
    \ingroup sugar
  */
  template<typename mtype,
	   typename diploid_t = std::pair<typename multiloc_glist_t<mtype>::iterator,
					  typename multiloc_glist_t<mtype>::iterator> >
  using multiloc = sugar::multiloc<mtype,
				   multiloc_mlist_t<mtype>,
				   multiloc_glist_t<mtype>,
				   boost::container::vector<boost::container::vector<diploid_t>>,
				   boost::container::vector<mtype>,
				   boost::container::vector<unsigned>,
				   boost::unordered_set<double,boost::hash<double>,KTfwd::equal_eps> >;

  /*!
    \brief Single population, multilocus simulation with serialization.  Can be copied, etc.
    \ingroup sugar
  */
  template<typename mtype,
	   typename mwriter,
	   typename mreader,
	   typename diploid_t = std::pair<typename multiloc_glist_t<mtype>::iterator,
					  typename multiloc_glist_t<mtype>::iterator>,
	   typename diploid_writer_t = diploidIOplaceholder,
	   typename diploid_reader_t = diploidIOplaceholder>
  using multiloc_serialized = sugar::multiloc_serialized<mtype,mwriter,mreader,
							 multiloc_mlist_t<mtype>,
							 multiloc_glist_t<mtype>,
							 boost::container::vector<boost::container::vector<diploid_t>>,
							 boost::container::vector<mtype>,
							 boost::container::vector<unsigned>,
							 boost::unordered_set<double,boost::hash<double>,KTfwd::equal_eps>,
							 diploid_writer_t,
							 diploid_reader_t>;
}

#else

#include <list>
#include <vector>
#include <unordered_set>
#include <fwdpp/sugar/multiloc/multiloc.hpp>
namespace KTfwd
{
  template<typename mtype> using multiloc_mlist_t = std::list<mtype>;
  template<typename mtype> using multiloc_gamete_t = gamete_base<mtype,multiloc_mlist_t<mtype>>;
  template<typename mtype> using multiloc_glist_t = std::list<multiloc_gamete_t<mtype>>;

  /*!
    \brief Single population, multilocus simulation without serialization.  Cannot be copied, etc.
    See @ref md_md_sugar for rationale, etc.
    \ingroup sugar
  */
  template<typename mtype,
	   typename diploid_t = std::pair<typename multiloc_glist_t<mtype>::iterator,
					  typename multiloc_glist_t<mtype>::iterator> >
  using multiloc = sugar::multiloc<mtype,
				   multiloc_mlist_t<mtype>,
				   multiloc_glist_t<mtype>,
				   std::vector<std::vector<diploid_t>>,
				   std::vector<mtype>,
				   std::vector<unsigned>,
				   std::unordered_set<double,std::hash<double>,KTfwd::equal_eps>>;
  
  /*!
    \brief Single population, multilocus simulation with serialization.  Can be copied, etc.
    See @ref md_md_sugar for rationale, etc.
    \ingroup sugar
  */
  template<typename mtype,
	   typename mwriter,
	   typename mreader,
	   typename diploid_t = std::pair<typename multiloc_glist_t<mtype>::iterator,
					  typename multiloc_glist_t<mtype>::iterator>,
	   typename diploid_writer_t = diploidIOplaceholder,
	   typename diploid_reader_t = diploidIOplaceholder>
  using multiloc_serialized = sugar::multiloc_serialized<mtype,mwriter,mreader,
							 multiloc_mlist_t<mtype>,
							 multiloc_glist_t<mtype>,
							 std::vector<std::vector<diploid_t>>,
							 std::vector<mtype>,
							 std::vector<unsigned>,
							 std::unordered_set<double,std::hash<double>,KTfwd::equal_eps>,
							 diploid_writer_t,
							 diploid_reader_t>;
}
#endif

#endif
