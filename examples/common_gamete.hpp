#ifndef __FWDPP_EXAMPLES_COMMON_GAMETE_HPP__
#define __FWDPP_EXAMPLES_COMMON_GAMETE_HPP__

#include <iostream>

#if defined(HAVE_BOOST_VECTOR) && defined(HAVE_BOOST_LIST) && defined(HAVE_BOOST_UNORDERED_SET) && defined(HAVE_BOOST_POOL_ALLOC) && defined(HAVE_BOOST_HASH) && !defined(USE_STANDARD_CONTAINERS)
#include <boost/unordered_set.hpp> 
#include <boost/container/list.hpp> 
#include <boost/container/vector.hpp>
#include <boost/pool/pool_alloc.hpp>
#include <boost/functional/hash.hpp>
typedef boost::pool_allocator<mtype> mut_allocator;
typedef boost::container::list<mtype,mut_allocator > mlist;
typedef KTfwd::gamete_base<mtype,mlist> gtype;
typedef boost::pool_allocator<gtype> gam_allocator;
typedef boost::container::vector<gtype,gam_allocator > gvector;
typedef boost::unordered_set<double,boost::hash<double>,KTfwd::equal_eps > lookup_table_type;
/*
  We wish to do mutations under the infinitely-many sites assumption.  That means that
  a new mutation cannot appear at any previously-mutated site.  Here, we cheat a little
  and do not allow mutations at any sites that are currently polymorphic.

  We accomplish this via a lookup table of the current mutation positions.  The function object
  KTfwd::equal_eps is used as a replacement for std::operator==(double,double) in order to ensure
  that values differing by <= DBL_EPSILON (~10^-17 on most systems) are not allowed, as they cause problems.
 */
typedef boost::unordered_set<double,boost::hash<double>,KTfwd::equal_eps > lookup_table_type;
#else
#include <unordered_set>
#include <vector>
#include <list>
typedef std::list<mtype> mlist;
typedef KTfwd::gamete_base<mtype,mlist> gtype;
typedef std::vector<gtype> gvector;
/*
  We wish to do mutations under the infinitely-many sites assumption.  That means that
  a new mutation cannot appear at any previously-mutated site.  Here, we cheat a little
  and do not allow mutations at any sites that are currently polymorphic.

  We accomplish this via a lookup table of the current mutation positions.  The function object
  KTfwd::equal_eps is used as a replacement for std::operator==(double,double) in order to ensure
  that values differing by <= DBL_EPSILON (~10^-17 on most systems) are not allowed, as they cause problems.
 */
typedef std::unordered_set<double,std::hash<double>,KTfwd::equal_eps > lookup_table_type;
#endif

#ifdef USE_STANDARD_CONTAINERS

#else

#endif

#endif
