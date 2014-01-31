//  -*- C++ -*- 
#ifndef _DIPLOID_FUNCTIONS_TCC_
#define _DIPLOID_FUNCTIONS_TCC_

#include <utility>
#include <cassert>
#include <algorithm>
#include <iterator>
#include <functional>

#include <cmath>
#include <numeric>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/is_base_and_derived.hpp>
#include <boost/container/vector.hpp>

#include <fwdpp/util.hpp>
#include <fwdpp/insertion_policies.hpp>
#include <fwdpp/fwd_functional.hpp>

//Implementation details in the following files
#include <fwdpp/diploid_functions_gamete_based.tcc>
#include <fwdpp/diploid_functions_ind_based.tcc>
#include <fwdpp/diploid_functions_recombine_gametes.tcc>
  
#endif /* _DIPLOID_FUNCTIONS_TCC_ */
