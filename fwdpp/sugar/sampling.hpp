#ifndef __FWDPP_SUGAR_SAMPLING_HPP__
#define __FWDPP_SUGAR_SAMPLING_HPP__

#include <fwdpp/sampling_functions.hpp>
#include <fwdpp/sugar/singlepop.hpp>
#include <fwdpp/sugar/metapop.hpp>
#include <fwdpp/sugar/multiloc.hpp>

namespace KTfwd
{
  using sample_t = std::vector< std::pair<double,std::string> >;
  using sep_sample_t = std::pair<sample_t,sample_t>;

  //Single-region, single-deme
  template<typename poptype>
  sep_sample_t sample_sep_details( gsl_rng * r,
				   const poptype & p,
				   const unsigned nsam,
				   const bool removeFixed,
				   std::true_type)
  {
    return ms_sample_separate(r,&p.diploids,nsam,removeFixed);
  }

  //Multi-locus, single-deme
  template<typename poptype>
  sep_sample_t sample_sep_details( gsl_rng * r,
				   const poptype & p,
				   const unsigned nsam,
				   const bool removeFixed,
				   std::false_type)
  {
    return ms_sample_separate(r,&p.diploids,nsam,removeFixed);
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
    return ms_sample_separate(r,&p.diploids[deme],nsam,removeFixed);
  }



}

#endif
