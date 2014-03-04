//  -*- C++ -*- 
#ifndef __DIPLOID_FUNCTIONS_REC_GAMS_TCC__
#define __DIPLOID_FUNCTIONS_REC_GAMS_TCC__

namespace KTfwd
{
 //recombination for individual-based simulation
  template< typename iterator_type,
	    typename gamete_insertion_policy,
	    typename recombination_map,
	    typename list_type_allocator,
	    template<typename,typename> class list_type>
  unsigned recombine_gametes( gsl_rng * r,
			      const double & littler,
			      list_type< typename iterator_type::value_type,list_type_allocator > * gametes,
			      iterator_type & g1,
			      iterator_type & g2,
			      const recombination_map & mf,
			      const gamete_insertion_policy & gpolicy)
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
	for(unsigned i = 0 ; i < nbreaks ; ++i)
	  {
	    pos.push_back(mf());
	  }
	std::sort(pos.begin(),pos.end());
	pos.push_back(std::numeric_limits<double>::max());
	
	typename iterator_type::value_type new_gamete1(0u,
						       typename iterator_type::value_type::mutation_container(),
						       typename iterator_type::value_type::mutation_container()),
	  new_gamete2(new_gamete1);
	new_gamete1.mutations.reserve(g1->mutations.size()+g2->mutations.size());
	new_gamete1.smutations.reserve(g1->smutations.size()+g2->smutations.size());
	new_gamete2.mutations.reserve(g1->mutations.size()+g2->mutations.size());
	new_gamete2.smutations.reserve(g1->smutations.size()+g2->smutations.size());
	short SWITCH_I = 0,SWITCH_J = 0;
	size_t dummy = 0;
	typename iterator_type::value_type::mcont_iterator itr = g1->mutations.begin(),
	  jtr = g2->mutations.begin(),
	  itr_s = g1->smutations.begin(),
	  jtr_s = g2->smutations.begin();
	//pointer arithmetic over a range of pointers.  apologies...
	//typename gamete_type::mutation_container::iterator itr2;
	for( ; dummy < pos.size(); ++dummy)
	  {
	    //iterate over neutral mutations from parent i
	    for( ; itr < g1->mutations.end() && (*itr)->pos < pos[dummy] ;++itr)
	      {
		switch(SWITCH_I)
		  {
		  case 0:
#ifndef NDEBUG
		    for(unsigned t = 0 ; t < new_gamete1.mutations.size() ; ++t )
		      {
			assert(! mutation_at_pos()(*new_gamete1.mutations[t],(*itr)->pos) );
		      }
#endif
		    new_gamete1.mutations.push_back(*itr);
		    break;
		  case 1:
#ifndef NDEBUG
		    for(unsigned t = 0 ; t < new_gamete2.mutations.size() ; ++t )
		      {
			assert(! mutation_at_pos()(*new_gamete2.mutations[t],(*itr)->pos) );
		      }
#endif
		    new_gamete2.mutations.push_back(*itr);
		    break;
		  }
	      }
	    
	    //iterate over selected mutations from parent i
	    for( ; itr_s < g1->smutations.end() && (*itr_s)->pos < pos[dummy] ;++itr_s)
	      {
		switch(SWITCH_I)
		  {
		  case 0:
#ifndef NDEBUG
		    for(unsigned t = 0 ; t < new_gamete1.smutations.size() ; ++t )
		      {
			assert(! mutation_at_pos()(*new_gamete1.smutations[t],(*itr_s)->pos) );
		      }
#endif
		    new_gamete1.smutations.push_back(*itr_s);
		    break;
		  case 1:
#ifndef NDEBUG
		    for(unsigned t = 0 ; t < new_gamete2.smutations.size() ; ++t )
		      {
			assert(! mutation_at_pos()(*new_gamete2.smutations[t],(*itr_s)->pos) );
		      }
#endif
		    new_gamete2.smutations.push_back(*itr_s);
		    break;
		  }
	      }
	    
	    //iterate over neutral mutations from parent j
	    for( ; jtr < g2->mutations.end() && (*jtr)->pos < pos[dummy] ;++jtr)
	      {
		switch(SWITCH_J)
		  {
		  case 1:
#ifndef NDEBUG
		  for(unsigned t = 0 ; t < new_gamete1.mutations.size() ; ++t )
		    {
		      assert(! mutation_at_pos()(*new_gamete1.mutations[t],(*jtr)->pos) );
		    }
#endif
		  new_gamete1.mutations.push_back(*jtr);
		  break;
		case 0:
#ifndef NDEBUG
		  for(unsigned t = 0 ; t < new_gamete2.mutations.size() ; ++t )
		    {
		      assert(! mutation_at_pos()(*new_gamete2.mutations[t],(*jtr)->pos) );
		    }
#endif
		  new_gamete2.mutations.push_back(*jtr);
		  break;
		}
	    }
	  
	  //iterate over selected mutations from parent j
	  for( ; jtr_s < g2->smutations.end() && (*jtr_s)->pos < pos[dummy] ;++jtr_s)
	    {
	      switch(SWITCH_J)
		{
		case 1:
#ifndef NDEBUG
		  for(unsigned t = 0 ; t < new_gamete1.smutations.size() ; ++t )
		    {
		      assert(! mutation_at_pos()(*new_gamete1.smutations[t],(*jtr_s)->pos) );
		    }
#endif
		  new_gamete1.smutations.push_back(*jtr_s);
		  break;
		case 0:
#ifndef NDEBUG
		  for(unsigned t = 0 ; t < new_gamete2.smutations.size() ; ++t )
		    {
		      assert(! mutation_at_pos()(*new_gamete2.smutations[t],(*jtr_s)->pos) );
		    }
#endif
		  new_gamete2.smutations.push_back(*jtr_s);
		  break;
		}
	    }
	  SWITCH_I=!SWITCH_I;
	  SWITCH_J=!SWITCH_J;
	  }
	std::sort(new_gamete1.mutations.begin(),new_gamete1.mutations.end(),fake_less());
	std::sort(new_gamete1.smutations.begin(),new_gamete1.smutations.end(),fake_less());
	std::sort(new_gamete2.mutations.begin(),new_gamete2.mutations.end(),fake_less());
	std::sort(new_gamete2.smutations.begin(),new_gamete2.smutations.end(),fake_less());
	g1 = gpolicy(new_gamete1,gametes);
	g2 = gpolicy(new_gamete2,gametes);
      }
    return nbreaks;
  }
}

#endif
