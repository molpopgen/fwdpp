#ifndef __FWDPP_SUGAR_METAPOP_HPP__
#define __FWDPP_SUGAR_METAPOP_HPP__

#include <fwdpp/sugar/metapop/metapop.hpp>
#include <fwdpp/fwd_functional.hpp>
#include <list>
#include <vector>
#include <unordered_set>

namespace KTfwd
{
  template<typename mtype> using metapop_mlist_t = std::list<mtype>;
  template<typename mtype> using metapop_gamete_t = gamete_base<mtype,metapop_mlist_t<mtype>>;
  template<typename mtype> using metapop_glist_t = std::list<metapop_gamete_t<mtype> >;

  /*!
    \brief Single locus metapopulation simulation without serialization.  Cannot be copied, etc.
    See @ref md_md_sugar for rationale, etc.
    \ingroup sugar
  */
  template<typename mtype,
	   typename diploid_t = std::pair<typename metapop_glist_t<mtype>::iterator,
					  typename metapop_glist_t<mtype>::iterator> >
  using metapop = sugar::metapop<mtype,
				 metapop_mlist_t<mtype>,
				 metapop_glist_t<mtype>,
				 std::vector<diploid_t>,
				 std::vector<std::vector<diploid_t> >,
				 std::vector<mtype>,
				 std::vector<uint_t>,
				 std::unordered_set<double,std::hash<double>,KTfwd::equal_eps>>;

  /*!
    \brief Single locus metapopulation simulation with serialization.  Can be copied, etc.
    See @ref md_md_sugar for rationale, etc.
    \ingroup sugar
  */
  template<typename mtype,
	   typename mwriter_t,
	   typename mreader_t,
	   typename diploid_t = std::pair<typename metapop_glist_t<mtype>::iterator,
					  typename metapop_glist_t<mtype>::iterator>,
	   typename diploid_writer_t = diploidIOplaceholder,
	   typename diploid_reader_t = diploidIOplaceholder>
  using metapop_serialized = sugar::metapop_serialized<mtype,mwriter_t,mreader_t,
						       metapop_mlist_t<mtype>,
						       metapop_glist_t<mtype>,
						       std::vector<diploid_t>,
						       std::vector<std::vector<diploid_t> >,
						       std::vector<mtype>,
						       std::vector<uint_t>,
						       std::unordered_set<double,std::hash<double>,KTfwd::equal_eps>,
						       diploid_writer_t,
						       diploid_reader_t>;
}
#endif
