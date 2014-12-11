//  -*- C++ -*- 
#ifndef __DIPLOID_FUNCTIONS_REC_GAMS_TCC__
#define __DIPLOID_FUNCTIONS_REC_GAMS_TCC__

#include <fwdpp/internal/recombination_common.hpp>

namespace KTfwd
{
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
			      const recombination_map & mf)
  {
    assert( g1 != gametes->end() );
    assert( g2 != gametes->end() );
    
    typedef typename iterator_type::value_type gtype;
    typedef typename gtype::mutation_container gtype_mcont;

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
	gtype new_gamete1(0u,gtype_mcont(),gtype_mcont()),
	  new_gamete2(new_gamete1);
	
 	new_gamete1.mutations.reserve(g1->mutations.size()+g2->mutations.size());
 	new_gamete1.smutations.reserve(g1->smutations.size()+g2->smutations.size());
 	new_gamete2.mutations.reserve(g1->mutations.size()+g2->mutations.size());
 	new_gamete2.smutations.reserve(g1->smutations.size()+g2->smutations.size());
	
	fwdpp_internal::recombine_gametes(pos,g1,g2,new_gamete1,new_gamete2);
	
	auto current_end = gametes->end();
	bool f1 = false, f2 = false;
	for( auto itr = gametes->begin() ;
	     (!f1||!f2)&&itr != current_end ; ++itr )
	  {
	    if(!f1&&*itr == new_gamete1)
	      {
		g1=itr;
		f1=true;
	      }
	    if(!f2&&*itr == new_gamete2)
	      {
		g2=itr;
		f2=true;
	      }
	  }
	if(!f1)
	  {
	    g1=gametes->insert(gametes->end(),std::move(new_gamete1));
	  }
	if(!f2)
	  {
	    g2=gametes->insert(gametes->end(),std::move(new_gamete2));
	  }
      }
    return nbreaks;
  }
}

#endif
