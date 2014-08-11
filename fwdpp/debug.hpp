#ifndef __KTFWD_DEBUG_HPP__
#define __KTFWD_DEBUG_HPP__

#ifndef NDEBUG

#include <fwdpp/forward_types.hpp>
#include <type_traits>

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
    static_assert( std::is_base_of<gamete_base_type,gamete_type>::value ||
                   std::is_same<gamete_base_type,gamete_type>::value,
                   "gamete_type must be, or inherit from, KTfwd::gamete_base<mutation_type,mutation_list_type>" );
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
    static_assert( std::is_base_of<gamete_base_type,gamete_type>::value ||
                   std::is_same<gamete_base_type,gamete_type>::value,
                   "gamete_type must be, or inherit from, KTfwd::gamete_base<mutation_type,mutation_list_type>" );
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
