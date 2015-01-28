//  -*- C++ -*- 
#ifndef __FWDPP_SAMPLING_FUNCTIONS_TCC__
#define __FWDPP_SAMPLING_FUNCTIONS_TCC__

#include <fwdpp/fwd_functional.hpp>
#include <fwdpp/internal/ms_sampling.hpp>
#include <limits>
#include <iostream>

namespace KTfwd
{
  template< typename gamete_type,
	    typename vector_type_allocator,
	    template<typename,typename> class vector_type>
  std::vector<unsigned> sample(gsl_rng * r,
			       const vector_type<gamete_type,vector_type_allocator > & gametes,
			       const unsigned & n, const unsigned & N)
  {
    std::vector<double> freqs(gametes.size(),0);
    std::vector<unsigned> counts(gametes.size(),0);
    for(unsigned i=0;i<gametes.size();++i)
      {
	double f = double(gametes[i].n)/double(N);
	freqs[i]=f;
      }
    gsl_ran_multinomial(r,gametes.size(),n,&freqs[0],&counts[0]);
    return counts;
  }

  template< typename gamete_type,
	    typename vector_type_allocator,
	    template<typename,typename> class vector_type>
  std::vector<unsigned> sample_sfs(gsl_rng * r, 
				   const vector_type<gamete_type,vector_type_allocator > & gametes,
				   const unsigned & n, const unsigned & N)
  {
    std::vector<unsigned> counts = sample(r,gametes,n,N);
    std::map<double,unsigned> samplemuts;
    std::map<double,unsigned>::iterator itr;
    for(unsigned i=0;i<gametes.size();++i)
      {
	if(counts[i]>0)
	  {
	    for(unsigned j=0;j<gametes[i].mutations.size();++j)
	      {
		itr = samplemuts.find(gametes[i].mutations[j]->pos);
		if( itr == samplemuts.end() )
		  {
		    samplemuts[gametes[i].mutations[j]->pos] = counts[i];
		  }
		else
		  {
		    itr->second += counts[i];
		  }
	      }
	    for(unsigned j=0;j<gametes[i].smutations.size();++j)
	      {
		itr = samplemuts.find(gametes[i].smutations[j]->pos);
		if( itr == samplemuts.end() )
		  {
		    samplemuts[gametes[i].smutations[j]->pos] = counts[i];
		  }
		else
		  {
		    itr->second += counts[i];
		  }
	      }
	  }
      }
    std::vector<unsigned> samplesfs(n,0);
    for(itr=samplemuts.begin();itr!=samplemuts.end();++itr)
      {
	samplesfs[itr->second-1]++;
      }
    return samplesfs;
  }

  template< typename gamete_type,
	    typename vector_type_allocator,
	    template<typename,typename> class vector_type>
  std::vector< std::pair<double, std::string> > 
  ms_sample(gsl_rng * r,
	    const vector_type<gamete_type,vector_type_allocator > & gametes,
	    const unsigned & n, const unsigned & N,
	    bool remove_fixed)
  {
    auto separate = ms_sample_separate(r,gametes,n,N,remove_fixed);
    std::move(separate.second.begin(),separate.second.end(),
    	      std::back_inserter(separate.first));
    std::sort(separate.first.begin(),
    	      separate.first.end(),
    	      [](std::pair<double,std::string> lhs,
    		 std::pair<double,std::string> rhs) { return lhs.first < rhs.first; });
    return separate.first;
  }

  template< typename gamete_type,
	    typename vector_type_allocator,
	    template<typename,typename> class vector_type>
  std::pair< std::vector< std::pair<double, std::string> > ,
	     std::vector< std::pair<double, std::string> > >
  ms_sample_separate(gsl_rng * r,
		     const vector_type<gamete_type,vector_type_allocator > & gametes,
		     const unsigned & n, const unsigned & N,
		     bool remove_fixed)
  {
    std::vector<unsigned> counts = sample(r,gametes,n,N);	
    std::vector< std::pair<double, std::string> > rvsel,rvneut;
    std::vector< std::pair<double, std::string> >::iterator itr;

    std::function<bool(const std::pair<double,std::string> &, const double &)> sitefinder = [](const std::pair<double,std::string> & site,
											       const double & d ) 
      {
	return std::fabs(site.first-d) <= std::numeric_limits<double>::epsilon();
      };

    unsigned individual = 0;
    for(unsigned i=0;i<counts.size();++i)
      {
	for(unsigned j=0;j<counts[i];++j,++individual)
	  {
	    fwdpp_internal::update_sample_block(rvneut,gametes[i].mutations,individual,n,sitefinder,0,1);
	    fwdpp_internal::update_sample_block(rvsel,gametes[i].smutations,individual,n,sitefinder,0,1);
	  }
      }
    assert(individual==n);
		
    if(remove_fixed)
      {
	rvneut.erase( std::remove_if(rvneut.begin(),rvneut.end(),[&n]( const std::pair<double,std::string> & site ) {
	      return unsigned(std::count(site.second.begin(),site.second.end(),'1')) == n; 
		} ),
	  rvneut.end() );
	rvsel.erase( std::remove_if(rvsel.begin(),rvsel.end(),[&n]( const std::pair<double,std::string> & site ) {
	      return unsigned(std::count(site.second.begin(),site.second.end(),'1')) == n; 
		} ),
	  rvsel.end() );
      }
    if(!rvneut.empty())
      {
	std::sort(rvneut.begin(),rvneut.end(),
		  [](std::pair<double,std::string> lhs,
		     std::pair<double,std::string> rhs) { return lhs.first < rhs.first; });
      }
    if(!rvsel.empty())
      {
	std::sort(rvsel.begin(),rvsel.end(),
		  [](std::pair<double,std::string> lhs,
		     std::pair<double,std::string> rhs) { return lhs.first < rhs.first; });
      }
    return std::make_pair(rvneut,rvsel);
  }

  //SAMPLERS FOR INDIVIDUAL-BASED SIMULATIONS

  template<typename iterator_type,
	   typename allocator,
	   template<typename,typename> class vector_type >
  std::vector< std::pair<double,std::string> >
  ms_sample( gsl_rng * r,
	     const vector_type< std::pair<iterator_type,iterator_type>, allocator > * diploids,
	     const unsigned & n,
	     const bool & remove_fixed)
  {
    auto separate = ms_sample_separate(r,diploids,n,remove_fixed);
    std::move( separate.second.begin(), separate.second.end(), std::back_inserter(separate.first) );
    std::sort(separate.first.begin(),separate.first.end(),
	      [](std::pair<double,std::string> lhs,
			   std::pair<double,std::string> rhs) { return lhs.first < rhs.first; });
    return separate.first;
  }

  template<typename iterator_type,
	   typename allocator,
	   template<typename,typename> class vector_type >
  std::pair<std::vector< std::pair<double,std::string> >,
	    std::vector< std::pair<double,std::string> > >
  ms_sample_separate( gsl_rng * r,
		      const vector_type< std::pair<iterator_type,iterator_type>, allocator > * diploids,
		      const unsigned & n,
		      const bool & remove_fixed)
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
    for( unsigned i = 0 ; i < n/2 ; ++i )
      {
	typename decltype(dptr)::difference_type ind = decltype(ind)(gsl_ran_flat(r,0.,double(diploids->size())));
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

  //Individual-based sims, multilocus algorithm
 

  template<typename iterator_type,
	   typename allocator,
	   typename outer_allocator,
	   template<typename,typename> class vector_type,
	   template<typename,typename> class outer_vector_type>
  std::vector< std::pair<std::vector< std::pair<double,std::string> >,std::vector< std::pair<double,std::string> > > >
  ms_sample_separate( gsl_rng * r,
		      const outer_vector_type< vector_type< std::pair<iterator_type,iterator_type>, allocator >, outer_allocator > * diploids,
		      const unsigned & n,
		      const bool & remove_fixed)
  {
    //typedef std::vector< std::vector< std::pair<double,std::string> > > rvtype;
    using rvtype = std::vector< std::pair<std::vector< std::pair<double,std::string> > ,
					  std::vector< std::pair<double,std::string> > > >;
    typedef std::vector< std::pair<double,std::string> >::iterator rv_inner_itr;
    typedef vector_type< std::pair<iterator_type,iterator_type>, allocator > genotype;
    typedef outer_vector_type< genotype, outer_allocator > dip_ctr;
    typedef typename iterator_type::value_type::mcont_iterator mut_itr;

    rvtype rv( diploids->size() );

    std::vector< typename dip_ctr::size_type > individuals;
    for( unsigned i = 0 ; i < n/2  ; ++i )
      {
	individuals.push_back( typename dip_ctr::size_type( gsl_ran_flat(r,0,diploids->size()) ) );
      }

    std::function<bool(const std::pair<double,std::string> &, const double &)> sitefinder = [](const std::pair<double,std::string> & site,
											       const double & d ) 
      {
	return std::fabs(site.first-d) <= std::numeric_limits<double>::epsilon();
      };

    //Go over each indidivual's mutations and update the return value
    typename dip_ctr::const_iterator dbegin = diploids->begin();
    for( unsigned ind = 0 ; ind < individuals.size() ; ++ind )
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

 template<typename iterator_type,
	   typename allocator,
	   typename outer_allocator,
	   template<typename,typename> class vector_type,
	   template<typename,typename> class outer_vector_type>
  std::vector< std::vector< std::pair<double,std::string> > >
  ms_sample( gsl_rng * r,
	     const outer_vector_type< vector_type< std::pair<iterator_type,iterator_type>, allocator >, outer_allocator > * diploids,
	     const unsigned & n,
	     const bool & remove_fixed)
  {
    auto separate = ms_sample_separate(r,diploids,n,remove_fixed);
    std::vector< std::vector< std::pair<double,std::string> > > rv;
    for( unsigned i = 0 ; i < separate.size() ; ++i )
      {
	std::move( separate[i].second.begin(), separate[i].second.end(),
		   std::back_inserter(separate[i].first) );
	std::sort(separate[i].firstbegin(),separate[i].first.end(),
		  [](std::pair<double,std::string> lhs,
		     std::pair<double,std::string> rhs) { return lhs.first < rhs.first; });
	rv.emplace_back(std::move(separate[i].first));
      }
    return rv;
    // typedef std::vector< std::vector< std::pair<double,std::string> > > rvtype;
    // typedef std::vector< std::pair<double,std::string> >::iterator rv_inner_itr;
    // typedef vector_type< std::pair<iterator_type,iterator_type>, allocator > genotype;
    // typedef outer_vector_type< genotype, outer_allocator > dip_ctr;
    // typedef typename iterator_type::value_type::mcont_iterator mut_itr;

    // rvtype rv( diploids->size() );

    // std::vector< typename dip_ctr::size_type > individuals;
    // for( unsigned i = 0 ; i < n/2  ; ++i )
    //   {
    // 	individuals.push_back( typename dip_ctr::size_type( gsl_ran_flat(r,0,diploids->size()) ) );
    //   }

    // std::function<bool(const std::pair<double,std::string> &, const double &)> sitefinder = [](const std::pair<double,std::string> & site,
    // 											       const double & d ) 
    //   {
    // 	return std::fabs(site.first-d) <= std::numeric_limits<double>::epsilon();
    //   };

    // //Go over each indidivual's mutations and update the return value
    // typename dip_ctr::const_iterator dbegin = diploids->begin();
    // for( unsigned ind = 0 ; ind < individuals.size() ; ++ind )
    //   {
    // 	unsigned rv_count=0;
    // 	for( typename genotype::const_iterator locus = (dbegin+ind)->begin() ; 
    // 	     locus < (dbegin+ind)->end() ; ++locus, ++rv_count )
    // 	  {
    // 	    //finally, we can go over mutations
    // 	    fwdpp_internal::update_sample_block(rv[rv_count],locus->first->mutations,ind,n,sitefinder);
    // 	    fwdpp_internal::update_sample_block(rv[rv_count],locus->first->smutations,ind,n,sitefinder);
    // 	    fwdpp_internal::update_sample_block(rv[rv_count],locus->second->mutations,ind,n,sitefinder,1);
    // 	    fwdpp_internal::update_sample_block(rv[rv_count],locus->second->smutations,ind,n,sitefinder,1);
    // 	  }
    //   }
  
    // if( remove_fixed )
    //   {
    // 	for( unsigned i = 0 ; i < rv.size() ; ++i )
    // 	  {
    // 	    if( ! rv[i].empty() )
    // 	      {
    // 		rv[i].erase( std::remove_if(rv[i].begin(),rv[i].end(),[&n]( const std::pair<double,std::string> & site ) {
    // 		      return unsigned(std::count(site.second.begin(),site.second.end(),'1')) == n; 
    // 		    } ),
    // 		  rv[i].end() );
    // 	      }
    // 	  }
    //   }
    // //sort on position
    // for( unsigned i = 0 ; i < rv.size() ; ++i )
    //   {
    // 	std::sort(rv[i].begin(),rv[i].end(),
    // 		  [](std::pair<double,std::string> lhs,
    // 		     std::pair<double,std::string> rhs) { return lhs.first < rhs.first; });
    //   }
    // return rv;
  }
}

#endif
