//  -*- C++ -*- 
#ifndef __DIPLOID_FUNCTIONS_REC_GAMS_TCC__
#define __DIPLOID_FUNCTIONS_REC_GAMS_TCC__

#include <iostream>
#include <fwdpp/algorithm.hpp>
namespace KTfwd
{
  namespace
  {
    struct rec_gamete_updater
    {
      template<typename itr_type,
	       typename cont_type>
      inline bool operator()( itr_type & i, cont_type * m1, cont_type * m2,
			      const unsigned & SWITCH, const double & val ) const
      {
	if( i->pos < val )
	  {
	    if( SWITCH )
	      {
		m1->emplace_back(i);
	      }
	    else
	      {
		m2->emplace_back(i);
	      }
	    return true;
	  }
	return false;
      }
    };
    
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
			      const recombination_map & mf)
  {
    assert( g1 != gametes->end() );
    assert( g2 != gametes->end() );
    
    //Identify cases where recombination cannot result in changed gametes, and get out quick
    if(g1 == g2 ) return 0;
    unsigned nm1=g1->mutations.size()+g1->smutations.size();
    unsigned nm2=g2->mutations.size()+g2->smutations.size();
    if((std::min(nm1,nm2)==0 && std::max(nm1,nm2)==1)) return 0;
    
    unsigned nbreaks = (littler > 0) ? gsl_ran_poisson(r,littler) : 0.;
    
    if( nbreaks )
      {
#ifndef USE_STANDARD_CONTAINERS
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
	typename iterator_type::value_type new_gamete1(0u,
						       typename iterator_type::value_type::mutation_container(),
						       typename iterator_type::value_type::mutation_container()),
	  new_gamete2(new_gamete1);

	new_gamete1.mutations.reserve(g1->mutations.size()+g2->mutations.size());
	new_gamete1.smutations.reserve(g1->smutations.size()+g2->smutations.size());
	new_gamete2.mutations.reserve(g1->mutations.size()+g2->mutations.size());
	new_gamete2.smutations.reserve(g1->smutations.size()+g2->smutations.size());
	
	typename iterator_type::value_type::mcont_const_iterator itr = g1->mutations.cbegin(),
	  jtr = g2->mutations.cbegin(),
	  itr_s = g1->smutations.cbegin(),
	  jtr_s = g2->smutations.cbegin(),
	  itr_e = g1->mutations.cend(),
	  itr_s_e = g1->smutations.cend(),
	  jtr_e = g2->mutations.cend(),
	  jtr_s_e = g2->smutations.cend();
	short SWITCH = 0;
	rec_gamete_updater UPDATER;
	for(const auto dummy : pos)
	  {
	    for_each_if( itr, itr_e,
			 std::bind(UPDATER,std::placeholders::_1,
				   &new_gamete2.mutations,&new_gamete1.mutations,SWITCH,dummy));
	    for_each_if( itr_s, itr_s_e,
			 std::bind(UPDATER,std::placeholders::_1,
				   &new_gamete2.smutations,&new_gamete1.smutations,SWITCH,dummy));
	    for_each_if( jtr, jtr_e,
			 std::bind(UPDATER,std::placeholders::_1,
				   &new_gamete1.mutations,&new_gamete2.mutations,SWITCH,dummy));
	    for_each_if( jtr_s, jtr_s_e,
			 std::bind(UPDATER,std::placeholders::_1,
				   &new_gamete1.smutations,&new_gamete2.smutations,SWITCH,dummy));
 	    SWITCH=!SWITCH;
 	  }

	//Below is the code used through fwdpp 0.2.4
// 	for(unsigned dummy=0 ; dummy < pos.size(); ++dummy)
// 	  {
// 	    //iterate over neutral mutations from parent i
// 	    for( ; itr < g1->mutations.end() && (*itr)->pos < pos[dummy] ;++itr)
// 	      {
// 		switch(SWITCH)
// 		  {
// 		  case 0:
// #ifndef NDEBUG
// 		    for(unsigned t = 0 ; t < new_gamete1.mutations.size() ; ++t )
// 		      {
// 			assert(! mutation_at_pos()(*new_gamete1.mutations[t],(*itr)->pos) );
// 		      }
// #endif
// 		    new_gamete1.mutations.push_back(*itr);
// 		    break;
// 		  case 1:
// #ifndef NDEBUG
// 		    for(unsigned t = 0 ; t < new_gamete2.mutations.size() ; ++t )
// 		      {
// 			assert(! mutation_at_pos()(*new_gamete2.mutations[t],(*itr)->pos) );
// 		      }
// #endif
// 		    new_gamete2.mutations.push_back(*itr);
// 		    break;
// 		  }
// 	      }
	    
// 	    //iterate over selected mutations from parent i
// 	    for( ; itr_s < g1->smutations.end() && (*itr_s)->pos < pos[dummy] ;++itr_s)
// 	      {
// 		switch(SWITCH)
// 		  {
// 		  case 0:
// #ifndef NDEBUG
// 		    for(unsigned t = 0 ; t < new_gamete1.smutations.size() ; ++t )
// 		      {
// 			assert(! mutation_at_pos()(*new_gamete1.smutations[t],(*itr_s)->pos) );
// 		      }
// #endif
// 		    new_gamete1.smutations.push_back(*itr_s);
// 		    break;
// 		  case 1:
// #ifndef NDEBUG
// 		    for(unsigned t = 0 ; t < new_gamete2.smutations.size() ; ++t )
// 		      {
// 			assert(! mutation_at_pos()(*new_gamete2.smutations[t],(*itr_s)->pos) );
// 		      }
// #endif
// 		    new_gamete2.smutations.push_back(*itr_s);
// 		    break;
// 		  }
// 	      }
	    
// 	    //iterate over neutral mutations from parent j
// 	    for( ; jtr < g2->mutations.end() && (*jtr)->pos < pos[dummy] ;++jtr)
// 	      {
// 		switch(SWITCH)
// 		  {
// 		  case 1:
// #ifndef NDEBUG
// 		    for(unsigned t = 0 ; t < new_gamete1.mutations.size() ; ++t )
// 		      {
// 			assert(! mutation_at_pos()(*new_gamete1.mutations[t],(*jtr)->pos) );
// 		      }
// #endif
// 		    new_gamete1.mutations.push_back(*jtr);
// 		    break;
// 		  case 0:
// #ifndef NDEBUG
// 		    for(unsigned t = 0 ; t < new_gamete2.mutations.size() ; ++t )
// 		      {
// 			assert(! mutation_at_pos()(*new_gamete2.mutations[t],(*jtr)->pos) );
// 		      }
// #endif
// 		    new_gamete2.mutations.push_back(*jtr);
// 		    break;
// 		  }
// 	      }
	  
// 	    //iterate over selected mutations from parent j
// 	    for( ; jtr_s < g2->smutations.end() && (*jtr_s)->pos < pos[dummy] ;++jtr_s)
// 	      {
// 		switch(SWITCH)
// 		  {
// 		  case 1:
// #ifndef NDEBUG
// 		    for(unsigned t = 0 ; t < new_gamete1.smutations.size() ; ++t )
// 		      {
// 			assert(! mutation_at_pos()(*new_gamete1.smutations[t],(*jtr_s)->pos) );
// 		      }
// #endif
// 		    new_gamete1.smutations.push_back(*jtr_s);
// 		    break;
// 		  case 0:
// #ifndef NDEBUG
// 		    for(unsigned t = 0 ; t < new_gamete2.smutations.size() ; ++t )
// 		      {
// 			assert(! mutation_at_pos()(*new_gamete2.smutations[t],(*jtr_s)->pos) );
// 		      }
// #endif
// 		    new_gamete2.smutations.push_back(*jtr_s);
// 		    break;
// 		  }
// 	      }
// 	    SWITCH=!SWITCH;
// 	  }

	std::sort(new_gamete1.mutations.begin(),new_gamete1.mutations.end(),fake_less());
	std::sort(new_gamete1.smutations.begin(),new_gamete1.smutations.end(),fake_less());
	std::sort(new_gamete2.mutations.begin(),new_gamete2.mutations.end(),fake_less());
	std::sort(new_gamete2.smutations.begin(),new_gamete2.smutations.end(),fake_less());
	typename list_type< typename iterator_type::value_type,list_type_allocator >::iterator current_end = gametes->end();
	bool f1 = false, f2 = false;
	for( typename list_type< typename iterator_type::value_type,list_type_allocator >::iterator itr = gametes->begin() ;
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
