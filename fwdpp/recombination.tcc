//  -*- C++ -*- 
#ifndef __FWDPP_RECOMBINATION_TCC__
#define __FWDPP_RECOMBINATION_TCC__

#include <vector>
#include <fwdpp/internal/recombination_common.hpp>
#include <fwdpp/internal/recycling.hpp>
namespace KTfwd
{  
  template< typename iterator_type,
	    typename list_type_allocator,
	    typename glookup_t,
	    typename queue_t,
	    template<typename,typename> class list_type>
  unsigned recombine_gametes( std::vector<double> & pos,
			      list_type< typename iterator_type::value_type,list_type_allocator > * gametes,
			      iterator_type & g1,
			      iterator_type & g2,
			      glookup_t & gamete_lookup,
			      queue_t & gamete_recycling_bin,
			      typename iterator_type::value_type::mutation_container & neutral,
			      typename iterator_type::value_type::mutation_container & selected )
  {    
    assert( g1 != gametes->end() );
    assert( g2 != gametes->end() );
    assert( std::is_sorted(pos.begin(),pos.end()) );
    assert( *(pos.end()-1) == std::numeric_limits<double>::max() );

    //We defer clearing all the way to this point
    neutral.clear();
    selected.clear();
    fwdpp_internal::recombine_gametes(pos,g1,g2,neutral,selected);

    //typename iterator_type::value_type ng(0u,neutral,selected);

    //Lookup table method modified in 0.3.5.  Result is faster simulations with selection.
    //auto lookup = gamete_lookup.lookup(ng);
    auto lookup = gamete_lookup.lookup(neutral,selected);
    if( lookup.first != lookup.second ) 
      {
	//Then we have to search through lookup.second
	auto itr = std::find_if(lookup.first,
				lookup.second,
				[&neutral,&selected]( typename glookup_t::inner_t & __p) {
				  return (__p.second->mutations == neutral &&
				  __p.second->smutations == selected);
				});
       if( itr == lookup.second )
	 {
	   fwdpp_internal::recycle_gamete(g1,gametes,gamete_recycling_bin,gamete_lookup,neutral,selected);
	 } 
	else 
	{
	  g1 = itr->second;
	}
      } 
    else
      {
	/*
	  There is no gamete in gametes with the number of neutral AND selected mutations,
	  and therefore the gamete is novel
	*/
	fwdpp_internal::recycle_gamete(g1,gametes,gamete_recycling_bin,gamete_lookup,neutral,selected);
      }
    return unsigned(pos.size()-1);
  }

  //recombination for individual-based simulation
  template< typename iterator_type,
	    typename recombination_map,
	    typename list_type_allocator,
	    typename glookup_t,
	    typename queue_t,
	    template<typename,typename> class list_type>
  unsigned recombine_gametes( gsl_rng * r,
			      const double & littler,
			      list_type< typename iterator_type::value_type,list_type_allocator > * gametes,
			      iterator_type & g1,
			      iterator_type & g2,
			      glookup_t & gamete_lookup,
			      queue_t & gamete_recycling_bin,
			      typename iterator_type::value_type::mutation_container & neutral,
			      typename iterator_type::value_type::mutation_container & selected,
			      const recombination_map & mf)
  {
    assert( g1 != gametes->end() );
    assert( g2 != gametes->end() );
    
    //Identify cases where recombination cannot result in changed gametes, and get out quick
    if(g1 == g2 ) return 0;
    auto nm1=g1->mutations.size()+g1->smutations.size();
    auto nm2=g2->mutations.size()+g2->smutations.size();
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
	return recombine_gametes(pos,gametes,g1,g2,gamete_lookup,gamete_recycling_bin,neutral,selected);
      }
    return 0;
  }
}

#endif
