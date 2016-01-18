/*!
  \file demography.hpp

  \brief Low-level functions for demographic events
*/

#ifndef FWDPP_DEMOGRAPHY_HPP
#define FWDPP_DEMOGRAPHY_HPP

#include <cassert>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <set>
#include <gsl/gsl_rng.h>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/internal/sample_diploid_helpers.hpp>
#include <fwdpp/internal/demography_details.hpp>

namespace KTfwd
{
  template<typename mcont_t,
	   typename mcount_t,
	   typename gcont_t,
	   typename vdipvector_t>
  int copy_deme(const mcont_t & mutations,
		mcount_t & mcounts,
		gcont_t & gametes,
		vdipvector_t & diploids,
		const size_t i)
  {
    if(i>=diploids.size()) return -1;
    diploids.emplace_back(diploids[i]);
    //We've added a deme, so we need to go through it, and
    //update gamete counts accordingly
    for( auto & dip : diploids[diploids.size()-1] )
      {
	gametes[dip.first].n++;
	gametes[dip.second].n++;
      }
    fwdpp_internal::process_gametes(gametes,mutations,mcounts);
    return 0;
  }

  template<typename vdipvector_t>
  int merge_demes(vdipvector_t & diploids,
		  size_t i,
		  size_t j)
  {
    if( i>=diploids.size() || j >= diploids.size() ) return -1;
    if(i==j) return 1;
    if(i>j) std::swap(i,j);
    std::move(diploids[j].begin(),diploids[j].end(),
	      std::back_inserter(diploids[i]));
    diploids.erase(diploids.begin()+j);
    return 0;
  }

  template<typename mcont_t,
	   typename mcount_t,
	   typename gcont_t,
	   typename vdipvector_t>
  int remove_deme(const mcont_t & mutations,
		  mcount_t & mcounts,
		  gcont_t & gametes,
		  vdipvector_t & diploids,
		  const size_t i)
  {
    if( i >= diploids.size() ) return -1;
    for( auto & dip : diploids[i] ) //update gamete counts
      {
	gametes[dip.first].n--;
	gametes[dip.second].n--;
      }
    diploids.erase(diploids.begin()+i);
    KTfwd::fwdpp_internal::process_gametes(gametes,mutations,mcounts); //update mutation counts
    return 0;
  }

  template<typename vdipvector_t>
  int swap_demes(vdipvector_t & diploids,
		 const size_t i, const size_t j)
  {
    if(i>=diploids.size()||j>=diploids.size()) return -1;
    if(i==j) return 1;
    std::swap(diploids[i],diploids[j]);
    return 0;
  }

  template<typename mcont_t,
	   typename mcount_t,
	   typename gcont_t,
	   typename vdipvector_t>
  int split_deme( gsl_rng * r,
		  const mcont_t & mutations,
		  mcount_t & mcounts,
		  gcont_t & gametes,
		  vdipvector_t & diploids,
		  const size_t i,
		  const uint_t N_new,
		  const bool replacement = false)
  {
    if(i>=diploids.size()) return -1;
    if(N_new >= diploids[i].size()) return 1;
    return fwdpp_internal::split_deme_details(r,mutations,mcounts,gametes,diploids,i,N_new,replacement);
  }

  template<typename mcont_t,
	   typename mcount_t,
	   typename gcont_t,
	   typename vdipvector_t>
  int admix_demes(gsl_rng * r,
		  const mcont_t & mutations,
		  mcount_t & mcounts,
		  gcont_t & gametes,
		  vdipvector_t & diploids,
		  const size_t i,
		  const size_t j,
		  const double pi,
		  const uint_t N_new,
		  const bool replacement = false)
  {
    if(i>=diploids.size()||j>=diploids.size()) return -1;
    if(pi<0.||pi>=1.) return 1;

    uint_t N_from_i = std::round(pi*double(N_new)), N_from_j = N_new - N_from_i;

    if(!replacement) //check for logic errors in input
      {
	if(N_from_i >= diploids[i].size()) return 1;
	if(N_from_j >= diploids[j].size()) return 1;
      }
    
    //add a new deme
    diploids.emplace_back(typename vdipvector_t::value_type());
    //get reference to new deme
    auto & new_deme = diploids[diploids.size()-1];
    new_deme.reserve(N_new);
    const auto & parental_deme_i = diploids[i];

    //get list of individuals from deme i
    auto indlist = fwdpp_internal::sample_individuals(r,diploids[i].size(),N_from_i,replacement);

    //add individuals from deme i into new deme
    for(const auto & ind : indlist) new_deme.push_back(parental_deme_i[ind]);

    //get list of individuals from deme j
    indlist = fwdpp_internal::sample_individuals(r,diploids[j].size(),N_from_j,replacement);

    //update reference to parental deme
    const auto & parental_deme_j = diploids[j];
    
    //add individuals from deme i into new deme
    for(const auto & ind : indlist) new_deme.push_back(parental_deme_j[ind]);

    //update gamete counts due to new deme
    for(const auto & dip : new_deme)
      {
	gametes[dip.first].n++;
	gametes[dip.second].n++;
      }

    fwdpp_internal::process_gametes(gametes,mutations,mcounts);
    return 0;
  }
  
}

#endif
