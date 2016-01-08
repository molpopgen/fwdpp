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
   */
  template<typename T, typename cT>
  inline typename cT::iterator insert_at_end(  T && t,  cT * ct )
  {
    return ct->insert(ct->end(),std::forward<T>(t));
  }

  struct emplace_back
  {
    template<typename T,typename cT>
    inline std::size_t operator()(T && t,cT & ct) const
    {
      ct.emplace_back(std::forward<T>(t));
      return ct.size()-1;
    }
  };

  /*! \brief     An insertion policy 
    
    Useful when mutation models can give rise to gametes identical to those currently existing in the population. 
    This will happen under finite-sites models.
    
    If the new gamete (ng) is not found in gametes, then it is inserted at the end.
    
    Otherwise, if it is found, the count of that gamete in the population is inremented by ng.n (the count of new gamete).
    
  */
  template<typename T, typename cT>
  inline typename cT::iterator update_if_exists_insert( T && t,  cT * ct )
  {
    typename cT::iterator itr = std::find(ct->begin(),ct->end(),t);
    if(itr == ct->end())
      {
	return insert_at_end(std::forward<T>(t),ct);
      }
    else
      {
	itr->n += t.n;
      }
    return itr;
  }

  template<typename T, typename cT>
  inline typename cT::iterator insert_if_not_found( T && t,  cT * ct )
  {
    typename cT::iterator itr = std::find(ct->begin(),ct->end(),t);
    if(itr == ct->end())
      {
	return insert_at_end(std::forward<T>(t),ct);
      }
    return itr;
  }

  template<typename T, typename cT>
  inline typename cT::iterator insert_if_not_found( T && t,  cT * ct, const unsigned & n )
  {
    typename cT::iterator itr = std::find(ct->begin(),ct->end(),t);
    if(itr == ct->end())
      {
	auto __ii = insert_at_end(std::forward<T>(t),ct);
	__ii->n=n;
	return __ii;
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
  inline void update_if_exists_push(  T && t,  cT * ct )
  {
    typename cT::iterator itr = std::find(ct->begin(),ct->end(),t);
    if(itr == ct->end())
      {
	push_at_end(std::forward<T>(t),ct);
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
  inline typename cT::iterator insert_new_or_fail( T && t,  cT * ct )
  {
    typename cT::iterator itr = std::find(ct->begin(),ct->end(),t);
    if(itr == ct->end())
      {
	return insert_at_end(std::forward<T>(t),ct);
      }
    //else, return end() to indicate no insertion
    return ct->end();
  }
}
#endif /* _INSERTION_POLICIES_HPP_ */
