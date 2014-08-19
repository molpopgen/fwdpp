#ifndef __FWDPP_ALGORITHM_HPP__
#define __FWDPP_ALGORITHM_HPP__

namespace KTfwd
{
  template<typename itr,typename function>
  itr for_each_if( itr & i, const itr & j, const function & f)
  {
    //if(i==j) return i;
    bool stillgoing = true;
    //for( ; i != j && stillgoing ; )
    while( i!=j&&stillgoing )
      {
	stillgoing = f(*i);
	if(stillgoing) ++i; else return i;
      }
    return i;  
  }
}

#endif
