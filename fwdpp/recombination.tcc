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
	    assert(gametes[itr->second].mutations==neutral);
	    assert(gametes[itr->second].smutations==selected);
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

  template<typename gcont_t,
	   typename mcont_t,
	   typename lookup_t,
	   typename recbin_t,
	   typename recpol_t>
  std::size_t recombination(gcont_t & gametes,
			    lookup_t & gamete_lookup,
			    recbin_t & gamete_recycling_bin,
			    typename gcont_t::value_type::mutation_container & neutral,
			    typename gcont_t::value_type::mutation_container & selected,
			    const recpol_t & rec_pol,
			    const std::size_t g1,
			    const std::size_t g2,
			    const mcont_t & mutations)
  {
    if(g1==g2) return g1;
    auto nm1=gametes[g1].mutations.size()+gametes[g1].smutations.size();
    auto nm2=gametes[g2].mutations.size()+gametes[g2].smutations.size();
    if((std::min(nm1,nm2)==0 && std::max(nm1,nm2)==1)) return g1;
    auto pos = rec_pol(gametes[g1],gametes[g2],mutations);
    if(pos.empty()) return g1;
    return recombine_gametes(pos,
			     gametes,
			     mutations,g1,g2,gamete_lookup,
			     gamete_recycling_bin,
			     neutral,selected);
  }
  
  template<typename gcont_t,
	   typename mcont_t,
	   typename lookup_t,
	   typename recbin_t,
	   typename recpol_t>
  std::size_t recombination(gcont_t & gametes,
			    lookup_t & gamete_lookup,
			    recbin_t & gamete_recycling_bin,
			    typename gcont_t::value_type::mutation_container & neutral,
			    typename gcont_t::value_type::mutation_container & selected,
			    const recpol_t & rec_pol,
			    const std::size_t g1,
			    const std::size_t g2,
			    const mcont_t & mutations,
			    unsigned * nrec)
  {
    static_assert( traits::valid_rec_model<recpol_t,typename gcont_t::value_type,mcont_t>::value,
		   "type recpol_t is not a valid recombination policy" );
    if(g1==g2) {*nrec=0;return g1;}
    auto nm1=gametes[g1].mutations.size()+gametes[g1].smutations.size();
    auto nm2=gametes[g2].mutations.size()+gametes[g2].smutations.size();
    if((std::min(nm1,nm2)==0 && std::max(nm1,nm2)==1)) {*nrec=0;return g1;}
    auto pos = rec_pol(gametes[g1],gametes[g2],mutations);
    if(pos.empty()) {*nrec=0;return g1;};
    assert(pos.back()==std::numeric_limits<double>::max());
    *nrec = pos.size()-1;
    return recombine_gametes(pos,
			     gametes,
			     mutations,g1,g2,gamete_lookup,
			     gamete_recycling_bin,
			     neutral,selected);
  }
}

#endif
