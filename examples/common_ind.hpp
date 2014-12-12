#ifndef __FWDPP_EXAMPLES_COMMON_IND_HPP__
#define __FWDPP_EXAMPLES_COMMON_IND_HPP__

#include <config.h>
#include <iostream>

#if defined(HAVE_BOOST_VECTOR) && defined(HAVE_BOOST_LIST) && defined(HAVE_BOOST_UNORDERED_SET) && defined(HAVE_BOOST_POOL_ALLOC) && defined(HAVE_BOOST_HASH) && !defined(USE_STANDARD_CONTAINERS)
#include <boost/container/vector.hpp>
#include <boost/container/list.hpp>
#include <boost/pool/pool_alloc.hpp>
#include <boost/unordered_set.hpp>
#include <boost/functional/hash.hpp>
typedef boost::pool_allocator<mtype> mut_allocator;
typedef boost::container::list<mtype,mut_allocator > mlist;
typedef KTfwd::gamete_base<mtype,mlist> gtype;
typedef boost::pool_allocator<gtype> gam_allocator;
typedef boost::container::list<gtype,gam_allocator > glist;
typedef boost::unordered_set<double,boost::hash<double>,KTfwd::equal_eps > lookup_table_type;
#else
#include <unordered_set>
#include <vector>
#include <list>
typedef std::list<mtype> mlist;
typedef KTfwd::gamete_base<mtype,mlist> gtype;
typedef std::list<gtype> glist;
typedef std::unordered_set<double,std::hash<double>,KTfwd::equal_eps > lookup_table_type;
#endif

#endif
