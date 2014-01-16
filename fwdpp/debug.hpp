#ifndef __KTFWD_DEBUG_HPP__
#define __KTFWD_DEBUG_HPP__

#ifndef NDEBUG

#include <fwdpp/forward_types.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/is_base_and_derived.hpp>

namespace KTfwd
{

  /*! \brief Returns true if the sum of counts in gametes equals twoN, false otherwise
    Returns true if the sum of counts in gametes equals twoN, false otherwise
   */
  template<typename vector_type_allocator,
	   typename gamete_type,
	   template<typename,typename> class vector_type>
  bool check_sum(const vector_type<gamete_type,vector_type_allocator> & gametes, const unsigned & twoN)
  {
    typedef gamete_base< typename gamete_type::mutation_type, typename gamete_type::mutation_list_type > gamete_base_type;
    BOOST_STATIC_ASSERT( (boost::is_base_and_derived<gamete_base_type,
			  gamete_type>::value) || (boost::is_same<gamete_base_type,gamete_type >::value) );
    unsigned check=0;
    for(typename vector_type<gamete_type,vector_type_allocator>::const_iterator i=gametes.begin();i!=gametes.end();++i)
      {
	check+=i->n;
      }
    return (check == twoN);
  }

  /*! \brief Returns true if the sum of counts in gametes equals twoN, false otherwise
    Returns true if the sum of counts in gametes equals twoN, false otherwise
   */
  template<typename vector_type_allocator,
	   typename gamete_type,
	   template<typename,typename> class vector_type>
  bool check_sum(const vector_type<gamete_type,vector_type_allocator> * gametes, const unsigned & twoN)
  {
    typedef gamete_base< typename gamete_type::mutation_type, typename gamete_type::mutation_list_type > gamete_base_type;
    BOOST_STATIC_ASSERT( (boost::is_base_and_derived<gamete_base_type,
			  gamete_type>::value) || (boost::is_same<gamete_base_type,gamete_type >::value) );
    unsigned check=0;
    for(typename vector_type<gamete_type,vector_type_allocator>::const_iterator i=gametes->begin();i!=gametes->end();++i)
      {
	check+=i->n;
      }
    return (check == twoN);
  }

}

#endif
#endif
