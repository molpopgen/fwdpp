//  -*- C++ -*- 
#ifndef __DIPLOID_FUNCTIONS_REC_GAMS_TCC__
#define __DIPLOID_FUNCTIONS_REC_GAMS_TCC__

#include <fwdpp/internal/recombination_common.hpp>

namespace KTfwd
{
  template< typename iterator_type,
	    typename list_type_allocator,
	    typename vector_type_allocator,
	    template<typename,typename> class vector_type,
	    template<typename,typename> class list_type>
  unsigned recombine_gametes( const vector_type< double, vector_type_allocator > & pos,
			      list_type< typename iterator_type::value_type,list_type_allocator > * gametes,
			      iterator_type & g1,
			      iterator_type & g2,
			      typename iterator_type::value_type::mutation_container & neutral,
			      typename iterator_type::value_type::mutation_container & selected )
  {
    assert( g1 != gametes->end() );
    assert( g2 != gametes->end() );
    assert( std::is_sorted(pos.begin(),pos.end()) );
    assert( *(pos.end()-1) == std::numeric_limits<double>::max() );

    fwdpp_internal::recombine_gametes(pos,g1,g2,neutral,selected);

    typename iterator_type::value_type ng(0u,neutral,selected);
    //IDEA: FUTURE: seems better to do a backwards search, right?
    auto itr = std::find(gametes->begin(),gametes->end(),ng);
    if(itr != gametes->end() )
      {
	g1 = itr;
      }
    else
      {
	g1 = gametes->emplace(gametes->end(),std::move(ng));
      }
    return pos.size()-1;
  }

  //recombination for individual-based simulation
  template< typename iterator_type,
	    typename recombination_map,
	    typename list_type_allocator,
	    template<typename,typename> class list_type>
  unsigned recombine_gametes( gsl_rng * r,
			      const double & littler,
			      list_type< typename iterator_type::value_type,list_type_allocator > * gametes,
			      iterator_type & g1,
			      iterator_type & g2,
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
#if defined(HAVE_BOOST_VECTOR) && !defined(USE_STANDARD_CONTAINERS)
        boost::container::vector<double> pos;
#else
        std::vector<double> pos;
#endif
	pos.reserve(nbreaks+1);
	for(unsigned i = 0 ; i < nbreaks ; ++i)
	  {
	    pos.emplace_back(mf());
	  }
	std::sort(pos.begin(),pos.end());
	pos.emplace_back(std::numeric_limits<double>::max());
	return recombine_gametes(pos,gametes,g1,g2,neutral,selected);
      }
    return 0;
  }
}

#endif
