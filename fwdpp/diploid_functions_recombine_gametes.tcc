//  -*- C++ -*- 
#ifndef __DIPLOID_FUNCTIONS_REC_GAMS_TCC__
#define __DIPLOID_FUNCTIONS_REC_GAMS_TCC__

#include <fwdpp/internal/recombination_common.hpp>

namespace KTfwd
{
  template< typename iterator_type,
	    typename list_type_allocator,
	    typename vector_type_allocator,
	    typename glookup_t,
	    template<typename,typename> class vector_type,
	    template<typename,typename> class list_type>
  unsigned recombine_gametes( const vector_type< double, vector_type_allocator > & pos,
			      list_type< typename iterator_type::value_type,list_type_allocator > * gametes,
			      iterator_type & g1,
			      iterator_type & g2,
			      glookup_t & gamete_lookup,
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

    typename iterator_type::value_type ng(0u,neutral,selected);
    //0.3.3
#ifndef VECTOR_GLOOKUP
    auto pitr = gamete_lookup.equal_range(neutral.size()+selected.size());
    auto itr = std::find_if(pitr.first,pitr.second,[&ng](const typename glookup_t::value_type & __p) {
	return *__p.second == ng;
      });
        if(itr == pitr.second)
      {
	g1 = gametes->emplace(gametes->end(),std::move(ng));
	gamete_lookup.insert(std::make_pair(g1->mutations.size()+g1->smutations.size(),g1));
      }
    else
      {
	g1 = itr->second;
      }
#else
    typename std::pair<std::int32_t,iterator_type> dummy(neutral.size()+selected.size(),g1);
    auto pitr = std::equal_range(gamete_lookup.begin(),gamete_lookup.end(),dummy,
     				 [](const typename std::pair<std::int32_t,iterator_type> & __p1, const typename std::pair<std::int32_t,iterator_type> & __p2 )
     				 {
     				   return __p1.first < __p2.first;
     				 });
    auto itr = std::find_if(pitr.first,pitr.second,[&ng](const typename glookup_t::value_type & __p) {
	return *__p.second == ng;
      });
    if(itr == pitr.second)
      {
	g1 = gametes->emplace(gametes->end(),std::move(ng));
	//By def'n, pitr.second is the place to insert...
	gamete_lookup.insert(pitr.second,std::make_pair(g1->mutations.size()+g1->smutations.size(),g1));
      }
    else
      {
	g1 = itr->second;
      }
#endif
    // std::cerr << std::distance(pitr.first,pitr.second) << ' ' << gamete_lookup.size() << ' '
    // 	      << (std::find_if(pitr.first,pitr.second,[&ng](typename glookup_t::value_type & __p) {
    // 		    return *__p.second == ng;
    // 		  }) == pitr.second)
    // 	      << '\n';
    /*
    //NOTE: backwards searches have little effect on speed, apparently...
    auto itr=std::find(gametes->begin(),gametes->end(),ng);
    if(itr!=gametes->end())g1=itr;
    else g1 = gametes->emplace(gametes->end(),std::move(ng));
    */
    return pos.size()-1;
  }

  //recombination for individual-based simulation
  template< typename iterator_type,
	    typename recombination_map,
	    typename list_type_allocator,
	    typename glookup_t,
	    template<typename,typename> class list_type>
  unsigned recombine_gametes( gsl_rng * r,
			      const double & littler,
			      list_type< typename iterator_type::value_type,list_type_allocator > * gametes,
			      iterator_type & g1,
			      iterator_type & g2,
			      glookup_t & gamete_lookup,
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
	return recombine_gametes(pos,gametes,g1,g2,gamete_lookup,neutral,selected);
      }
    return 0;
  }
}

#endif
