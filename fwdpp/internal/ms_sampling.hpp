#ifndef __FWDPP_INTERNAL_MS_SAMPLING_HPP__
#define __FWDPP_INTERNAL_MS_SAMPLING_HPP__

namespace KTfwd
{
  namespace fwdpp_internal
  {
    inline void remove_no_derived( sample_t * block )
    {
      block->erase( std::remove_if(block->begin(),block->end(),
				   [](std::pair<double,std::string> & p) {
				     return unsigned(std::count(p.second.begin(),p.second.end(),'0')) == p.second.size();
				   }), block->end() );
    }
    //Used when nsam is odd.  We just clip off the last individual
    inline void trim_last( sample_t * block )
    {
      std::for_each( block->begin(),block->end(),
		     []( std::pair<double,std::string> & p ) {
		       if(!p.second.empty())
			 {
			   //remove last character
			   p.second.erase(p.second.end()-1);
			 }
		     } );
      remove_no_derived(block);
    }

    template< typename mcont_t,
	      typename pos_finder>
    void update_sample_block(sample_t & block,
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
    sep_sample_t
    ms_sample_separate_single_deme( const dipvector_t * diploids,
				    const std::vector<unsigned> & diplist,
				    const unsigned n,
				    const bool & remove_fixed )
    {
      sep_sample_t rv;
      sample_t::iterator itr;

      std::function<bool(const sample_site_t &, const double &)> sitefinder = [](const sample_site_t & site,
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
	  fwdpp_internal::update_sample_block( rv.first,(dptr+ind)->first->mutations,i,2*diplist.size(),sitefinder);
	  fwdpp_internal::update_sample_block( rv.first,(dptr+ind)->second->mutations,i,2*diplist.size(),sitefinder,1);
	  fwdpp_internal::update_sample_block( rv.second,(dptr+ind)->first->smutations,i,2*diplist.size(),sitefinder);
	  fwdpp_internal::update_sample_block( rv.second,(dptr+ind)->second->smutations,i,2*diplist.size(),sitefinder,1);
	}
      if(remove_fixed&&!rv.first.empty())
	{
	  rv.first.erase( std::remove_if(rv.first.begin(),rv.first.end(),[&diplist]( const sample_site_t & site ) {
		return unsigned(std::count(site.second.begin(),site.second.end(),'1')) == diplist.size();
	      } ),
	    rv.first.end() );
	}
      if(!rv.first.empty())
	{
	  std::sort(rv.first.begin(),rv.first.end(),
		    [](const sample_site_t & lhs,
		       const sample_site_t & rhs) { return lhs.first < rhs.first; });
	}
      if(remove_fixed&&!rv.second.empty())
	{
	  rv.second.erase( std::remove_if(rv.second.begin(),rv.second.end(),[&diplist]( const std::pair<double,std::string> & site ) {
		return unsigned(std::count(site.second.begin(),site.second.end(),'1')) == 2*diplist.size();
	      } ),
	    rv.second.end() );
	}
      if(!rv.second.empty())
	{
	  std::sort(rv.second.begin(),rv.second.end(),
		    [](const sample_site_t & lhs,
		       const sample_site_t & rhs) { return lhs.first < rhs.first; });
	}
      //Deal w/odd sample sizes
      if(n%2 != 0.)
	{
	  trim_last(&rv.first);
	  trim_last(&rv.second);
	}
      return rv;
    }


    template<typename dipvector_t>
    std::vector<sep_sample_t >
    ms_sample_separate_mlocus( const dipvector_t * diploids,
			       const std::vector<unsigned> & diplist,
			       const unsigned & n,
			       const bool & remove_fixed)
    {
      using rvtype = std::vector<sep_sample_t>;
      using genotype = typename dipvector_t::value_type;

      rvtype rv( diploids->size() );

      std::function<bool(const sample_site_t &, const double &)> sitefinder = [](const sample_site_t & site,
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
	      fwdpp_internal::update_sample_block(rv[rv_count].first,locus->first->mutations,ind,2*diplist.size(),sitefinder);
	      fwdpp_internal::update_sample_block(rv[rv_count].second,locus->first->smutations,ind,2*diplist.size(),sitefinder);
	      fwdpp_internal::update_sample_block(rv[rv_count].first,locus->second->mutations,ind,2*diplist.size(),sitefinder,1);
	      fwdpp_internal::update_sample_block(rv[rv_count].second,locus->second->smutations,ind,2*diplist.size(),sitefinder,1);
	    }
	}

      if( remove_fixed )
	{
	  for( unsigned i = 0 ; i < rv.size() ; ++i )
	    {
	      rv[i].first.erase( std::remove_if(rv[i].first.begin(),rv[i].first.end(),[&diplist]( const sample_site_t & site ) {
		    return unsigned(std::count(site.second.begin(),site.second.end(),'1')) == 2*diplist.size();
		  } ),
		rv[i].first.end() );
	      rv[i].second.erase( std::remove_if(rv[i].second.begin(),rv[i].second.end(),[&diplist]( const sample_site_t & site ) {
		    return unsigned(std::count(site.second.begin(),site.second.end(),'1')) == 2*diplist.size();
		  } ),
		rv[i].second.end() );
	    }
	}
      //sort on position
      for( unsigned i = 0 ; i < rv.size() ; ++i )
	{
	  std::sort(rv[i].first.begin(),rv[i].first.end(),
		    [](const sample_site_t & lhs,
		       const sample_site_t & rhs) { return lhs.first < rhs.first; });
	  std::sort(rv[i].second.begin(),rv[i].second.end(),
		    [](const sample_site_t & lhs,
		       const sample_site_t & rhs) { return lhs.first < rhs.first; });
	  //Deal w/odd sample sizes
	  if(n%2!=0.)
	    {

	      trim_last(&rv[i].first);
	      trim_last(&rv[i].second);
	    }
	}

      return rv;
    }
  }
}
#endif
