#ifndef __FWDPP_SUGAR_SAMPLING_HPP__
#define __FWDPP_SUGAR_SAMPLING_HPP__

#include <algorithm>
#include <stdexcept>
#include <fwdpp/type_traits.hpp>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/sampling_functions.hpp>
#include <fwdpp/sugar/singlepop.hpp>
#include <fwdpp/sugar/metapop.hpp>
#include <fwdpp/sugar/multiloc.hpp>
#include <fwdpp/sugar/sampling/sampling_details.hpp>

namespace KTfwd
{
  /*!
    Take a random sample of size 'nsam' from a population
   */
  template<typename poptype>
  sample_t sample( const gsl_rng * r,
		   const poptype & p,
		   const unsigned nsam,
		   const bool removeFixed)
  /*!
    Take a random sample of nsam chromosomes from a population

    \param r A random-number generator
    \param p A population
    \param nsam The sample size
    \param removeFixed Whether or not to remove variants present in all nsam chromosomes

    \return A vector of both neutral and non-neutral variants
  */
  {
    static_assert( (std::is_same<typename poptype::popmodel_t,sugar::SINGLEPOP_TAG>::value ||
		    std::is_same<typename poptype::popmodel_t,sugar::MULTILOCPOP_TAG>::value ),
		   "poptype must be SINGLEPOP_TAG or MULTILOCPOP_TAG"
		   );
    return sample_details(r,p,nsam,removeFixed,typename std::is_same<typename poptype::popmodel_t,sugar::SINGLEPOP_TAG>::type());
  }

  template<typename poptype>
  sep_sample_t sample_separate( const gsl_rng * r,
				const poptype & p,
				const unsigned nsam,
				const bool removeFixed)
  /*!
    Take a random sample of nsam chromosomes from a population

    \param r A random-number generator
    \param p A population
    \param nsam The sample size
    \param removeFixed Whether or not to remove variants present in all nsam chromosomes

    \return A pair of vectors.  The first element contains neutral variants.  The second contains non-neutral variants.
  */
  {
    static_assert( (std::is_same<typename poptype::popmodel_t,sugar::SINGLEPOP_TAG>::value ||
		    std::is_same<typename poptype::popmodel_t,sugar::MULTILOCPOP_TAG>::value ),
		   "poptype must be SINGLEPOP_TAG or MULTILOCPOP_TAG"
		   );
    return sample_sep_details(r,p,nsam,removeFixed,typename std::is_same<typename poptype::popmodel_t,sugar::SINGLEPOP_TAG>::type());
  }

  template<typename poptype>
  sample_t sample(const poptype & p,
		  const std::vector<unsigned> & individuals,
		  const bool removeFixed)
  /*!
    Take a non-random sample of diploids from a population

    \param p A population
    \param individuals The indexes of the diploids to sample
    \param removeFixed Whether or not to remove variants present in all nsam chromosomes

    \return A vector of both neutral and non-neutral variants
  */
  {
    static_assert( (std::is_same<typename poptype::popmodel_t,sugar::SINGLEPOP_TAG>::value ||
		    std::is_same<typename poptype::popmodel_t,sugar::MULTILOCPOP_TAG>::value ),
		   "poptype must be SINGLEPOP_TAG or MULTILOCPOP_TAG"
		   );
    if (individuals.empty())return sample_t();
    if( std::find_if(individuals.begin(),individuals.end(),[&p](const unsigned & u) {
	  return u >= p.diploids.size();
	}) != individuals.end() )
      {
	throw std::out_of_range("KTfwd::sample_separate: individual index out of range");
      }
    return sample_details(p,individuals,removeFixed,typename std::is_same<typename poptype::popmodel_t,sugar::SINGLEPOP_TAG>::type()); 
  }
  
  template<typename poptype>
  sep_sample_t sample_separate(const poptype & p,
			       const std::vector<unsigned> & individuals,
			       const bool removeFixed)
  /*!
    Take a non-random sample of nsam chromosomes from a population
    
    \param p A population
    \param individuals The indexes of the diploids to 'view' from the population
    \param nsam The sample size
    \param removeFixed Whether or not to remove variants present in all nsam chromosomes
    
    \return A pair of vectors.  The first element contains neutral variants.  The second contains non-neutral variants.
  */
  {
    static_assert( (std::is_same<typename poptype::popmodel_t,sugar::SINGLEPOP_TAG>::value ||
		    std::is_same<typename poptype::popmodel_t,sugar::MULTILOCPOP_TAG>::value ),
		   "poptype must be SINGLEPOP_TAG or MULTILOCPOP_TAG"
		   );
    if (individuals.empty())return sep_sample_t();
    if( std::find_if(individuals.begin(),individuals.end(),[&p](const unsigned & u) {
	  return u >= p.diploids.size();
	}) != individuals.end() )
      {
	throw std::out_of_range("KTfwd::sample_separate: individual index out of range");
      }
    return sample_sep_details(p,individuals,removeFixed,typename std::is_same<typename poptype::popmodel_t,sugar::SINGLEPOP_TAG>::type()); 
  }

  template<typename poptype>
  sample_t sample( const gsl_rng * r,
		   const poptype & p,
		   const unsigned deme,
		   const unsigned nsam,
		   const bool removeFixed )
  /*!
    Take a random sample of nsam chromosomes from a meta-population

    \param r A random-number generator
    \param p A population
    \param p the index of the deme to sample
    \param nsam The sample size
    \param removeFixed Whether or not to remove variants present in all nsam chromosomes

    \return A vector of both neutral and non-neutral variants
  */
  {
    static_assert( std::is_same<typename poptype::popmodel_t,sugar::METAPOP_TAG>::value,
		   "METAPOP_TAG required");
    if(deme >= p.diploids.size())
      {
	throw std::out_of_range("KTfwd::sample_separate: deme index out of range");
      }
    auto temp = ms_sample_separate(r,p.mutations,p.gametes,p.diploids[deme],nsam,removeFixed);
    auto rv = std::move(temp.first);
    std::move(temp.second.begin(),temp.second.end(),std::back_inserter(rv));
    finish_sample(rv,p.fixations,nsam,removeFixed,sugar::treat_neutral::ALL);
    /*
    std::sort(rv.begin(),rv.end(),[](const sample_site_t & a,
				     const sample_site_t & b){
		return a.first<b.first;
	      });
    if(!removeFixed)
      add_fixations(rv,p.fixations,nsam,sugar::treat_neutral::ALL);
    */
    return rv;
  }
  
  template<typename poptype>
  sep_sample_t sample_separate( const gsl_rng * r,
				const poptype & p,
				const unsigned deme,
				const unsigned nsam,
				const bool removeFixed )
  /*!
    Take a random sample of nsam chromosomes from a meta-population

    \param r A random-number generator
    \param p A population
    \param deme the index of the deme to sample
    \param nsam The sample size
    \param removeFixed Whether or not to remove variants present in all nsam chromosomes

    \return A pair of vectors.  The first element contains neutral variants.  The second contains non-neutral variants.
  */
  {
    static_assert( std::is_same<typename poptype::popmodel_t,sugar::METAPOP_TAG>::value,
		   "METAPOP_TAG required");
    if(deme >= p.diploids.size())
      {
	throw std::out_of_range("KTfwd::sample_separate: deme index out of range");
      }
    auto x = ms_sample_separate(r,p.mutations,p.gametes,p.diploids[deme],nsam,removeFixed);
    finish_sample(x.first,p.fixations,nsam,removeFixed,sugar::treat_neutral::NEUTRAL);
    finish_sample(x.second,p.fixations,nsam,removeFixed,sugar::treat_neutral::SELECTED);
    return x;
  }

  template<typename poptype>
  sample_t sample(const poptype & p,
		  const unsigned deme,
		  const std::vector<unsigned> & individuals,
		  const bool removeFixed )
  /*!
    Take a non-random sample of nsam chromosomes from a meta-population
    
    \param p A population
    \param p the index of the deme to sample
    \param nsam The sample size
    \param removeFixed Whether or not to remove variants present in all nsam chromosomes
    
    \return A vector of both neutral and non-neutral variants
  */
  {
    static_assert( std::is_same<typename poptype::popmodel_t,sugar::METAPOP_TAG>::value,
		   "METAPOP_TAG required");
    if(deme >= p.diploids.size())
      {
	throw std::out_of_range("KTfwd::sample_separate: deme index out of range");
      }
    if(individuals.empty()) return sample_t();
    for( const auto i : individuals )
      {
	if(i>=p.diploids[deme].size())
	  {
	    throw std::out_of_range("KTfwd::sample_separate: individual index out of range");
	  }
      }
    auto temp = fwdpp_internal::ms_sample_separate_single_deme(p.mutations,p.gametes,p.diploids[deme],individuals,2*individuals.size(),removeFixed);
    auto rv = std::move(temp.first);
    std::move(temp.second.begin(),temp.second.end(),std::back_inserter(rv));
    finish_sample(rv,p.fixations,2*individuals.size(),removeFixed,sugar::treat_neutral::ALL);
    /*
    std::sort(rv.begin(),rv.end(),[](const sample_site_t & a,
				     const sample_site_t & b){
		return a.first<b.first;
	      });
    if(!removeFixed)
      add_fixations(rv,p.fixations,2*individuals.size(),sugar::treat_neutral::ALL);
    */
    return rv;
  }
  
  template<typename poptype>
  sep_sample_t sample_separate(const poptype & p,
			       const unsigned deme,
			       const std::vector<unsigned> & individuals,
			       const bool removeFixed )
  /*!
    Take a non-random sample of nsam chromosomes from a meta-population

    \param p A population
    \param deme the index of the deme to sample
    \param nsam The sample size
    \param removeFixed Whether or not to remove variants present in all nsam chromosomes

    \return A pair of vectors.  The first element contains neutral variants.  The second contains non-neutral variants.
  */
  {
    static_assert( std::is_same<typename poptype::popmodel_t,sugar::METAPOP_TAG>::value,
		   "METAPOP_TAG required");
    if(deme >= p.diploids.size())
      {
	throw std::out_of_range("KTfwd::sample_separate: deme index out of range");
      }
    if(individuals.empty()) return sep_sample_t();
    for( const auto i : individuals )
      {
	if(i>=p.diploids[deme].size())
	  {
	    throw std::out_of_range("KTfwd::sample_separate: individual index out of range");
	  }
      }
    auto x = fwdpp_internal::ms_sample_separate_single_deme(p.mutations,p.gametes,p.diploids[deme],individuals,2*individuals.size(),removeFixed);
    finish_sample(x.first,p.fixations,2*individuals.size(),removeFixed,sugar::treat_neutral::NEUTRAL);
    finish_sample(x.second,p.fixations,2*individuals.size(),removeFixed,sugar::treat_neutral::SELECTED);
    return x;
  }
}

#endif
