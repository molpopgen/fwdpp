#ifndef __FWDPP_INTERNAL_MS_SAMPLING_HPP__
#define __FWDPP_INTERNAL_MS_SAMPLING_HPP__

namespace KTfwd
{
  namespace fwdpp_internal
  {
    template< typename mcont_t,
	      typename pos_finder>
    void update_sample_block(std::vector< std::pair<double,std::string> > & block,
			     const mcont_t & mutations,
			     const unsigned & i,
			     const unsigned & n,
			     const pos_finder & pf,
			     const unsigned & offset = 0,
			     const unsigned & scalar = 2)
    {
      for( auto mptr = mutations.cbegin() ; mptr != mutations.cend() ; ++mptr )
	{
	  double mutpos = (*mptr)->pos;
	  auto itr = std::find_if(block.begin(),block.end(),std::bind(pf,std::placeholders::_1,mutpos));
	  if( itr == block.end() )
	    {
	      block.push_back( std::make_pair(mutpos,std::string(n,'0')) );
	      block[block.size()-1].second[scalar*i+offset] = '1';
	    }
	  else
	    {
	      itr->second[scalar*i+offset]='1';
	    }
	}
    }

    template<typename dipvector_t>
    std::pair<std::vector< std::pair<double,std::string> >,
	      std::vector< std::pair<double,std::string> > >
    ms_sample_separate_single_deme( const dipvector_t * diploids,
				    const std::vector<unsigned> & diplist,
				    const unsigned n,
				    const bool & remove_fixed )
    {
      std::pair<std::vector< std::pair<double,std::string> >,
		std::vector< std::pair<double,std::string> > > rv;
      std::vector< std::pair<double, std::string> >::iterator itr;
  
      std::function<bool(const std::pair<double,std::string> &, const double &)> sitefinder = [](const std::pair<double,std::string> & site,
												 const double & d ) 
	{
	  return std::fabs(site.first-d) <= std::numeric_limits<double>::epsilon();
	};

      const auto dptr = diploids->cbegin();
      for( unsigned i = 0 ; i < diplist.size() ; ++i )
	{
	  unsigned ind = diplist[i];
	  assert(ind>=0);
	  assert( unsigned(ind) < diploids->size() );
	  fwdpp_internal::update_sample_block( rv.first,(dptr+ind)->first->mutations,i,n,sitefinder);
	  fwdpp_internal::update_sample_block( rv.first,(dptr+ind)->second->mutations,i,n,sitefinder,1);
	  fwdpp_internal::update_sample_block( rv.second,(dptr+ind)->first->smutations,i,n,sitefinder);
	  fwdpp_internal::update_sample_block( rv.second,(dptr+ind)->second->smutations,i,n,sitefinder,1);
	}
      if(remove_fixed&&!rv.first.empty())
	{
	  rv.first.erase( std::remove_if(rv.first.begin(),rv.first.end(),[&n]( const std::pair<double,std::string> & site ) {
		return unsigned(std::count(site.second.begin(),site.second.end(),'1')) == n; 
	      } ),
	    rv.first.end() );
	}
      if(!rv.first.empty())
	{
	  std::sort(rv.first.begin(),rv.first.end(),
		    [](std::pair<double,std::string> lhs,
		       std::pair<double,std::string> rhs) { return lhs.first < rhs.first; });
	}
      if(remove_fixed&&!rv.second.empty())
	{
	  rv.second.erase( std::remove_if(rv.second.begin(),rv.second.end(),[&n]( const std::pair<double,std::string> & site ) {
		return unsigned(std::count(site.second.begin(),site.second.end(),'1')) == n; 
	      } ),
	    rv.second.end() );
	}
      if(!rv.second.empty())
	{
	  std::sort(rv.second.begin(),rv.second.end(),
		    [](std::pair<double,std::string> lhs,
		       std::pair<double,std::string> rhs) { return lhs.first < rhs.first; });
	}
      return rv;
    }
  

    template<typename dipvector_t>
    std::vector<std::pair<std::vector< std::pair<double,std::string> >,
			  std::vector< std::pair<double,std::string> > > >
    ms_sample_separate_mlocus( const dipvector_t * diploids,
			       const std::vector<unsigned> & diplist,
			       const unsigned & n,
			       const bool & remove_fixed)
    {
      using rvtype = std::vector< std::pair<std::vector< std::pair<double,std::string> > ,
					    std::vector< std::pair<double,std::string> > > >;
      using genotype = typename dipvector_t::value_type;
      //using dipvector_t = outer_vector_type< genotype, outer_allocator >;

      rvtype rv( diploids->size() );

      std::function<bool(const std::pair<double,std::string> &, const double &)> sitefinder = [](const std::pair<double,std::string> & site,
												 const double & d ) 
	{
	  return std::fabs(site.first-d) <= std::numeric_limits<double>::epsilon();
	};

      //Go over each indidivual's mutations and update the return value
      typename dipvector_t::const_iterator dbegin = diploids->begin();
      for( unsigned ind = 0 ; ind < diplist.size() ; ++ind )
	{
	  unsigned rv_count=0;
	  for( typename genotype::const_iterator locus = (dbegin+ind)->begin() ; 
	       locus < (dbegin+ind)->end() ; ++locus, ++rv_count )
	    {
	      //finally, we can go over mutations
	      fwdpp_internal::update_sample_block(rv[rv_count].first,locus->first->mutations,ind,n,sitefinder);
	      fwdpp_internal::update_sample_block(rv[rv_count].second,locus->first->smutations,ind,n,sitefinder);
	      fwdpp_internal::update_sample_block(rv[rv_count].first,locus->second->mutations,ind,n,sitefinder,1);
	      fwdpp_internal::update_sample_block(rv[rv_count].second,locus->second->smutations,ind,n,sitefinder,1);
	    }
	}
  
      if( remove_fixed )
	{
	  for( unsigned i = 0 ; i < rv.size() ; ++i )
	    {
	      rv[i].first.erase( std::remove_if(rv[i].first.begin(),rv[i].first.end(),[&n]( const std::pair<double,std::string> & site ) {
		    return unsigned(std::count(site.second.begin(),site.second.end(),'1')) == n; 
		  } ),
		rv[i].first.end() );
	      rv[i].second.erase( std::remove_if(rv[i].second.begin(),rv[i].second.end(),[&n]( const std::pair<double,std::string> & site ) {
		    return unsigned(std::count(site.second.begin(),site.second.end(),'1')) == n; 
		  } ),
		rv[i].second.end() );
	    }
	}
      //sort on position
      for( unsigned i = 0 ; i < rv.size() ; ++i )
	{
	  std::sort(rv[i].first.begin(),rv[i].first.end(),
		    [](std::pair<double,std::string> lhs,
		       std::pair<double,std::string> rhs) { return lhs.first < rhs.first; });
	  std::sort(rv[i].second.begin(),rv[i].second.end(),
		    [](std::pair<double,std::string> lhs,
		       std::pair<double,std::string> rhs) { return lhs.first < rhs.first; });
	}
      return rv;
    }
  }
}
#endif
