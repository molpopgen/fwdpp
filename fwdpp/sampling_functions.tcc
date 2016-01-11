//  -*- C++ -*- 
#ifndef __FWDPP_SAMPLING_FUNCTIONS_TCC__
#define __FWDPP_SAMPLING_FUNCTIONS_TCC__

#include <fwdpp/fwd_functional.hpp>
#include <fwdpp/internal/ms_sampling.hpp>
#include <limits>
#include <type_traits>
#include <algorithm>

namespace KTfwd
{
  template< typename gamete_type,
	    typename allocator_t,
	    template<typename,typename> class container_type>
  std::vector<unsigned> sample(gsl_rng * r,
			       const container_type<gamete_type,allocator_t > & gametes,
			       const unsigned & n, const unsigned & N)
  {
    std::vector<double> freqs;
    std::vector<unsigned> counts(gametes.size(),0);
    std::for_each( gametes.begin(), gametes.end(), [&freqs,&N](const gamete_type & __g) {
	freqs.emplace_back( std::move( double(__g.n)/double(N) ) );
      } );
    gsl_ran_multinomial(r,gametes.size(),n,&freqs[0],&counts[0]);
    return counts;
  }

  template< typename gamete_type,
	    typename allocator_t,
	    template<typename,typename> class container_type>
  std::vector<unsigned> sample_sfs(gsl_rng * r, 
				   const container_type<gamete_type,allocator_t > & gametes,
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

  //SAMPLERS FOR INDIVIDUAL-BASED SIMULATIONS
  template<typename mcont_t,
	   typename gcont_t,
	   typename allocator,
	   typename diploid_geno_t,
	   template<typename,typename> class vector_type >
  typename std::enable_if< std::is_base_of<mutation_base,typename mcont_t::value_type>::value,
			   sample_t >::type
  ms_sample( gsl_rng * r,
	     const mcont_t & mutations,
	     const gcont_t & gametes,
	     const vector_type< diploid_geno_t, allocator > & diploids,
	     const unsigned & n,
	     const bool & remove_fixed)
  {
    auto separate = ms_sample_separate(r,mutations,gametes,diploids,n,remove_fixed);
    std::move( separate.second.begin(), separate.second.end(), std::back_inserter(separate.first) );
    std::sort(separate.first.begin(),separate.first.end(),
	      [](const sample_site_t & lhs,
		 const sample_site_t & rhs) { return lhs.first < rhs.first; });
    return separate.first;
  }

  template<typename mcont_t,
	   typename gcont_t,
	   typename allocator,
	   typename diploid_geno_t,
	   template<typename,typename> class vector_type >
  typename std::enable_if< std::is_base_of<mutation_base,typename mcont_t::value_type>::value,
			   sep_sample_t >::type
  ms_sample_separate( gsl_rng * r,
		      const mcont_t & mutations,
		      const gcont_t & gametes,
		      const vector_type< diploid_geno_t, allocator > & diploids,
		      const unsigned & n,
		      const bool & remove_fixed)
  {
    std::vector<unsigned> diplist;
    unsigned isodd = (n%2 != 0.) ? 1u : 0u;
    for( unsigned i = 0 ; i < n/2+isodd ; ++i )
      {
	diplist.push_back(std::vector<unsigned>::value_type(gsl_ran_flat(r,0.,double(diploids.size()))));
      }
    return fwdpp_internal::ms_sample_separate_single_deme(mutations,gametes,diploids,diplist,n,remove_fixed);
  }

  //Individual-based sims, multilocus algorithm
  template<typename mcont_t,
	   typename gcont_t,
	   typename diploid_geno_t,
	   typename allocator,
	   typename outer_allocator,
	   template<typename,typename> class vector_type,
	   template<typename,typename> class outer_vector_type>
  typename std::enable_if< std::is_base_of<mutation_base,typename mcont_t::value_type>::value,
			   std::vector<sep_sample_t > >::type
  ms_sample_separate( gsl_rng * r,
		      const mcont_t & mutations,
		      const gcont_t & gametes,
		      const outer_vector_type< vector_type< diploid_geno_t, allocator >, outer_allocator > & diploids,
		      const unsigned & n,
		      const bool & remove_fixed)
  {
    std::vector<unsigned> diplist;
    unsigned isodd = (n%2 != 0.) ? 1u : 0u;
    for( unsigned i = 0 ; i < n/2+isodd ; ++i )
      {
	diplist.push_back(unsigned(gsl_ran_flat(r,0.,double(diploids.size()))));
      }
    return fwdpp_internal::ms_sample_separate_mlocus(mutations,gametes,diploids,diplist,n,remove_fixed);
  }

  template<typename mcont_t,
	   typename gcont_t,
	   typename diploid_geno_t,
	   typename allocator,
	   typename outer_allocator,
	   template<typename,typename> class vector_type,
	   template<typename,typename> class outer_vector_type>
  typename std::enable_if< std::is_base_of<mutation_base,typename mcont_t::value_type>::value,
			   std::vector< sample_t > >::type
  ms_sample( gsl_rng * r,
	     const mcont_t & mutations,
	     const gcont_t & gametes,
	     const outer_vector_type< vector_type< diploid_geno_t, allocator >, outer_allocator > & diploids,
	     const unsigned & n,
	     const bool & remove_fixed)
  {
    auto separate = ms_sample_separate(r,mutations,gametes,diploids,n,remove_fixed);
    std::vector<sample_t> rv;
    for( unsigned i = 0 ; i < separate.size() ; ++i )
      {
	std::move( separate[i].second.begin(), separate[i].second.end(),
		   std::back_inserter(separate[i].first) );
	std::sort(separate[i].first.begin(),separate[i].first.end(),
		  [](const sample_site_t & lhs,
		     const sample_site_t & rhs) { return lhs.first < rhs.first; });
	rv.emplace_back(std::move(separate[i].first));
      }
    return rv;
  }
}

#endif
