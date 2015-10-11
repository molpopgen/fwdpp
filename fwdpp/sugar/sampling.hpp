#ifndef __FWDPP_SUGAR_SAMPLING_HPP__
#define __FWDPP_SUGAR_SAMPLING_HPP__

#include <stdexcept>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/sampling_functions.hpp>
#include <fwdpp/sugar/singlepop.hpp>
#include <fwdpp/sugar/metapop.hpp>
#include <fwdpp/sugar/multiloc.hpp>

/*
  TODO:

  1.  fxns to sample specific diploids
  2.  Move implementation details into subfolder or .tcc file
*/
namespace KTfwd
{
  using sample_t = std::vector< std::pair<double,std::string> >;
  using sep_sample_t = std::pair<sample_t,sample_t>;

  enum class treat_neutral {ALL,NEUTRAL,SELECTED};

  template<typename vec_mutation_t>
  void add_fixations( sample_t * sample,
		      const vec_mutation_t & fixations,
		      const unsigned nsam,
		      const treat_neutral treat )
  {
    for( const auto & f : fixations)
      {
	if( treat == treat_neutral::ALL )
	  {
	    sample->emplace_back( std::make_pair(f.pos,std::string(nsam,'1')) );
	  }
	else if (treat == treat_neutral::NEUTRAL && f.neutral ) //only add neutral mutations
	  {
	    sample->emplace_back( std::make_pair(f.pos,std::string(nsam,'1')) );
	  }
	else if (treat == treat_neutral::SELECTED && !f.neutral ) //only add selected mutations
	  {
	    sample->emplace_back( std::make_pair(f.pos,std::string(nsam,'1')) );
	  }
      }
  }

  //Single-region, single-deme
  template<typename poptype>
  sample_t sample_details( gsl_rng * r,
			   const poptype & p,
			   const unsigned nsam,
			   const bool removeFixed,
			   std::true_type)
  {
    sample_t rv =  ms_sample(r,&p.diploids,nsam,removeFixed);
    if(!removeFixed)
      add_fixations(&rv,p.fixations,nsam,treat_neutral::ALL);
    return rv;
  }

  //Multi-locus, single-deme
  template<typename poptype>
  sample_t sample_details( gsl_rng * r,
				   const poptype & p,
			   const unsigned nsam,
			   const bool removeFixed,
			   std::false_type)
  {
    sample_t rv = ms_sample(r,&p.diploids,nsam,removeFixed);
    if(!removeFixed)
      add_fixations(&rv,p.fixations,nsam,treat_neutral::ALL);
    return rv;
  }

  //Single-region, single-deme
  template<typename poptype>
  sep_sample_t sample_sep_details( gsl_rng * r,
				   const poptype & p,
				   const unsigned nsam,
				   const bool removeFixed,
				   std::true_type)
  {
    sep_sample_t rv = ms_sample_separate(r,&p.diploids,nsam,removeFixed);
    if(! removeFixed)
      {
	add_fixations(&rv.first,p.fixations,nsam,treat_neutral::NEUTRAL);
	add_fixations(&rv.second,p.fixations,nsam,treat_neutral::SELECTED);
      }
    return rv;
  }

  //Multi-locus, single-deme
  template<typename poptype>
  sep_sample_t sample_sep_details( gsl_rng * r,
				   const poptype & p,
				   const unsigned nsam,
				   const bool removeFixed,
				   std::false_type)
  {
    sep_sample_t rv =  ms_sample_separate(r,&p.diploids,nsam,removeFixed);
    if(! removeFixed)
      {
	add_fixations(&rv.first,p.fixations,nsam,treat_neutral::NEUTRAL);
	add_fixations(&rv.second,p.fixations,nsam,treat_neutral::SELECTED);
      }
    return rv;
  }

  template<typename poptype>
  sample_t sample( gsl_rng * r,
		   const poptype & p,
		   const unsigned nsam,
		   const bool removeFixed)
  {
    static_assert( (std::is_same<typename poptype::popmodel_t,sugar::SINGLEPOP_TAG>::value ||
		    std::is_same<typename poptype::popmodel_t,sugar::MULTILOCPOP_TAG>::value ),
		   "poptype must be SINGLEPOP_TAG or MULTILOCPOP_TAG"
		   );
    return sample_details(r,p,nsam,removeFixed,typename std::is_same<typename poptype::popmodel_t,sugar::SINGLEPOP_TAG>::type());
  }

  template<typename poptype>
  sep_sample_t sample_separate( gsl_rng * r,
				const poptype & p,
				const unsigned nsam,
				const bool removeFixed)
  {
    static_assert( (std::is_same<typename poptype::popmodel_t,sugar::SINGLEPOP_TAG>::value ||
		    std::is_same<typename poptype::popmodel_t,sugar::MULTILOCPOP_TAG>::value ),
		   "poptype must be SINGLEPOP_TAG or MULTILOCPOP_TAG"
		   );
    return sample_sep_details(r,p,nsam,removeFixed,typename std::is_same<typename poptype::popmodel_t,sugar::SINGLEPOP_TAG>::type());
  }

  template<typename poptype>
  sep_sample_t sample_separate( gsl_rng * r,
				const poptype & p,
				const unsigned deme,
				const unsigned nsam,
				const bool removeFixed )
  {
    static_assert( std::is_same<typename poptype::popmodel_t,sugar::METAPOP_TAG>::value,
		   "METAPOP_TAG required");
    if(deme >= p.diploids.size())
      {
	throw std::out_of_range("KTfwd::sample_separate: deme index out of range");
      }
    return ms_sample_separate(r,&p.diploids[deme],nsam,removeFixed);
  }



}

#endif
