#ifndef _INSERTION_POLICIES_HPP_
#define _INSERTION_POLICIES_HPP_

#include <algorithm>
#include <fwdpp/fwd_functional.hpp> 
namespace KTfwd
{
  /*! \brief     An insertion policy

    A common use of this policy will be when a 
    mutation model always gives rise to a gamete that does not
    currently exist in the population, such as the infinitely-many sites
    model of mutations.  For such models, it is sufficient to simply append 
    the new gamete to the existing vector.

    Note: it is silly to pass ng to this function with ng.n = 0.

    \deprecated
   */
  template<typename T, typename cT>
  inline typename cT::iterator insert_at_end( const T & t,  cT * ct )
  {
    return ct->insert(ct->end(),t);
  }

  /*! \brief     An insertion policy
    Works via push_back instead of insert, so there is no return value.

    A common use of this policy will be when a 
    mutation model always gives rise to a gamete that does not
    currently exist in the population, such as the infinitely-many sites
    model of mutations.  For such models, it is sufficient to simply append 
    the new gamete to the existing vector.

    Note: it is silly to pass ng to this function with ng.n = 0.

    \deprecated
   */
  template<typename T, typename cT>
  inline void push_at_end( const T & t,  cT * ct )
  {
    ct->push_back(t);
  }

  /*! \brief     An insertion policy 
    
    Useful when mutation models can give rise to gametes identical to those currently existing in the population. 
    This will happen under finite-sites models.
    
    If the new gamete (ng) is not found in gametes, then it is inserted at the end.
    
    Otherwise, if it is found, the count of that gamete in the population is inremented by ng.n (the count of new gamete).
    
  */
  template<typename T, typename cT>
  inline typename cT::iterator update_if_exists_insert( const T & t,  cT * ct )
  {
    typename cT::iterator itr = std::find(ct->begin(),ct->end(),t);
    if(itr == ct->end())
      {
	return insert_at_end(t,ct);
      }
    else
      {
	itr->n += t.n;
      }
    return itr;
  }

  template<typename T, typename cT>
  inline typename cT::iterator insert_if_not_found( const T & t,  cT * ct )
  {
    typename cT::iterator itr = std::find(ct->begin(),ct->end(),t);
    if(itr == ct->end())
      {
	return insert_at_end(t,ct);
      }
    return itr;
  }

  /*! \brief    An insertion policy 
    
    Useful when mutation models can give rise to gametes identical to those currently existing in the population. 
    This will happen under finite-sites models.
    
    If the new gamete (ng) is not found in gametes, then it is inserted at the end.
    
    Otherwise, if it is found, the count of that gamete in the population is inremented by ng.n (the count of new gamete).
    
    Note: it is silly to pass ng to this function with ng.n = 0.
  */
  template<typename T, typename cT>
  inline void update_if_exists_push( const T & t,  cT * ct )
  {
    typename cT::iterator itr = std::find(ct->begin(),ct->end(),t);
    if(itr == ct->end())
      {
	push_at_end(t,ct);
      }
    else
      {
	itr->n += t.n;
      }
  }

  /*! \brief An insertion policy

    If t is not found in ct, then insert t at end of ct, else
    return ct->end()
  */
  template<typename T, typename cT>
  inline typename cT::iterator insert_new_or_fail( const T & t,  cT * ct )
  {
    typename cT::iterator itr = std::find(ct->begin(),ct->end(),t);
    if(itr == ct->end())
      {
	return insert_at_end(t,ct);
      }
    //else, return end() to indicate no insertion
    return ct->end();
  }

  template<typename T, typename cT>
  inline typename cT::iterator insert_before_pos(const T & t, cT  * ct)
  {
    typename cT::iterator p = std::find_if(ct->begin(),ct->end(),boost::bind(greater_pos_object(),_1,t.pos));
    return ct->insert(p,t);
  }

  template<typename T, typename cT>
  inline typename cT::iterator insert_before_pos(const T & t, const double & pos, cT  * ct)
  {
    typename cT::iterator p = std::find_if(ct->begin(),ct->end(),boost::bind(greater_pos_object(),_1,pos));
    return ct->insert(p,t);
  }

  struct push_back_gamete
  /*! 
    \brief An insertion policy for gametes.  

    A common use of this policy will be when a 
    mutation model always gives rise to a gamete that does not
    currently exist in the population, such as the infinitely-many sites
    model of mutations.  For such models, it is sufficient to simply append 
    the new gamete to the existing vector.

    Note: it is silly to pass ng to this function with ng.n = 0.

    \deprecated
   */
  {
    template< typename gamete_type,
	      typename vector_type_allocator,
	      template<typename,typename> class vector_type>//,
    inline void operator()( const gamete_type & ng, 
			    vector_type<gamete_type,vector_type_allocator > * gametes)const 
    {
      gametes->push_back(ng);
    }
  };

  struct insert_unique
  /*! 
    \brief An insertion policy for gametes.

    Useful when mutation models can give rise to gametes identical to those currently existing in the population. 
    This will happen under finite-sites models.

    If the new gamete (ng) is not found in gametes, then it is inserted at the end.

    Otherwise, if it is found, the count of that gamete in the population is inremented by ng.n (the count of new gamete).

    Note: it is silly to pass ng to this function with ng.n = 0.

    Note: Legacy/deprecated
  */
  {
    template<typename gamete_type,
	     typename vector_type_allocator,
	     template<typename,typename> class vector_type>
    inline void operator()(const gamete_type & ng, 
			   vector_type<gamete_type,vector_type_allocator > * gametes) const
    {
      typedef typename  vector_type<gamete_type,vector_type_allocator >::iterator vtype_iterator;
      vtype_iterator itr=std::find(gametes->begin(),gametes->end(),ng);
      if(itr == gametes->end())
	{
	  gametes->push_back(ng);
	}
      else
	{
	  itr->n += ng.n;
	}
    }
  };

  /*!
    \brief Insertion policy for a mutation.

    Simply appends it to end of list_type mutations

    \deprecated
   */
  template<typename mutation_type,typename list_type>
  inline typename list_type::iterator insert_mutation_at_end(const mutation_type & m,list_type * mutations)
  {
    return (mutations->insert(mutations->end(),m));
  }
  
  /*!
    \deprecated
  */
  template<typename mutation_type,
	 typename list_type>
  inline typename list_type::iterator
  insert_unique_or_fail(const mutation_type & m,list_type * mutations)
  {
    typename list_type::iterator itr = std::find(mutations->begin(),mutations->end(),m);
    if( itr != mutations->end() ) return (itr);
    return (mutations->insert(mutations->end(),m));
  }
}
#endif /* _INSERTION_POLICIES_HPP_ */
