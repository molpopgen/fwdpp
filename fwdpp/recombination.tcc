//  -*- C++ -*- 
#ifndef __FWDPP_RECOMBINATION_TCC__
#define __FWDPP_RECOMBINATION_TCC__

#include <vector>
#include <fwdpp/internal/recombination_common.hpp>
#include <fwdpp/internal/recycling.hpp>
namespace KTfwd
{  
  template< typename vec_t,
	    typename gcont_t,
	    typename mcont_t,
	    typename glookup_t,
	    typename queue_t>
  std::size_t recombine_gametes( const vec_t & pos,
				 gcont_t & gametes,
				 const mcont_t & mutations,
				 const std::size_t g1,
				 const std::size_t g2,
				 glookup_t & gamete_lookup,
				 queue_t & gamete_recycling_bin,
				 std::vector<std::size_t> & neutral,
				 std::vector<std::size_t> & selected )
  {
    assert(g1 < gametes.size());
    assert(g2 < gametes.size());
    assert( std::is_sorted(pos.begin(),pos.end()) );
    assert( *(pos.end()-1) == std::numeric_limits<double>::max() );

    //We defer clearing all the way to this point
    neutral.clear();
    selected.clear();
    fwdpp_internal::recombine_gametes(pos,g1,g2,gametes,mutations,neutral,selected);

    //Lookup table method modified in 0.3.5.  Result is faster simulations with selection.
    auto lookup = gamete_lookup.lookup(neutral,selected,mutations);
    if( lookup.first != lookup.second ) 
      {
	//Then we have to search through lookup.second
	auto itr = std::find_if(lookup.first,
				lookup.second,
				[&gametes,&neutral,&selected]( typename glookup_t::inner_t & __p) {
				  return (gametes[__p.second].mutations == neutral &&
					  gametes[__p.second].smutations == selected);
				});
       if( itr == lookup.second )
	 {
	   return fwdpp_internal::recycle_gamete(gametes,mutations,gamete_recycling_bin,gamete_lookup,neutral,selected);
	 } 
	else 
	{
	  return itr->second;
	}
      } 
    else
      {
	/*
	  There is no gamete in gametes with the number of neutral AND selected mutations,
	  and therefore the gamete is novel
	*/
	return fwdpp_internal::recycle_gamete(gametes,mutations,gamete_recycling_bin,gamete_lookup,neutral,selected);
      }
    return g1;
  }

  //recombination for individual-based simulation
  template< typename recombination_map,
	    typename gcont_t,
	    typename mcont_t,
	    typename glookup_t,
	    typename queue_t>
  std::size_t recombine_gametes( gsl_rng * r,
				 const double & littler,
				 gcont_t & gametes,
				 const mcont_t & mutations,
				 const std::size_t g1,
				 const std::size_t g2,
				 glookup_t & gamete_lookup,
				 queue_t & gamete_recycling_bin,
				 std::vector<std::size_t> & neutral,
				 std::vector<std::size_t> & selected,
				 const recombination_map & mf)
  {
    assert( g1 < gametes.size() );
    assert( g2 < gametes.size() );
    
    //Identify cases where recombination cannot result in changed gametes, and get out quick
    if(g1 == g2 ) return g1;
    auto nm1=gametes[g1].mutations.size()+gametes[g1].smutations.size();
    auto nm2=gametes[g2].mutations.size()+gametes[g2].smutations.size();
    if((std::min(nm1,nm2)==0 && std::max(nm1,nm2)==1)) return 0;
    
    unsigned nbreaks = (littler > 0) ? gsl_ran_poisson(r,littler) : 0u;
    
    if( nbreaks )
      {
        std::vector<double> pos;
	pos.reserve(nbreaks+1);
	for(unsigned i = 0 ; i < nbreaks ; ++i)
	  {
	    pos.emplace_back(mf());
	  }
	std::sort(pos.begin(),pos.end());
	pos.emplace_back(std::numeric_limits<double>::max());
	return recombine_gametes(pos,gametes,mutations,g1,g2,gamete_lookup,gamete_recycling_bin,neutral,selected);
      }
    return g1;
  }
}

#endif
