///! \file regions.hpp
#ifndef __FWDPP_EXTENSIONS_REGIONS_HPP__
#define __FWDPP_EXTENSIONS_REGIONS_HPP__

#include <fwdpp/sugar/infsites.hpp>
#include <fwdpp/internal/gsl_discrete.hpp>
#include <fwdpp/extensions/callbacks.hpp>

namespace KTfwd
{
  namespace extensions
  {
    struct discrete_mut_model
    /*!
      Class allowing the simulation of discrete variation
      in mutation models along a region.

      Implemented in terms of KTfwd::popgenmut.

      Intended to be used with the callbacks in fwdpp/extensions/callbacks.hpp
     */
    {
      using result_type = KTfwd::popgenmut;
      std::vector<double> nbeg,nend,sbeg,send;
      std::vector< shmodel > shmodels;
      KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr nlookup,slookup;

      /*!
	\param __nbeg Positions of beginnings of 'neutral' regions
	\param __nend Positions of ends of 'neutral' regions
	\param nweights Weights on 'neutra'l regions
	\param __sbeg Positions of beginnings of 'selected' regions
	\param __send Positions of ends of 'selected' regions
	\param sweights Weights on 'selected' regions
	\param __shmodels Vector of KTfwd::experimenta::shmodel
       */
      discrete_mut_model(std::vector<double> __nbeg,
			 std::vector<double> __nend,
			 std::vector<double>  nweights, //the weights
			 std::vector<double> __sbeg,
			 std::vector<double> __send,
			 const std::vector<double> & sweights, //the weights
			 std::vector<shmodel> __shmodels) : nbeg(std::move(__nbeg)),nend(std::move(__nend)),
							    sbeg(std::move(__sbeg)),send(std::move(__send)),
							    shmodels(std::move(__shmodels))
      {
	if(nweights.size())
	  nlookup = KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(gsl_ran_discrete_preproc(nweights.size(),&nweights[0]));
	if(sweights.size())
	  slookup = KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(gsl_ran_discrete_preproc(sweights.size(),&sweights[0]));
      }

      /*!
	Helper function.  Not to be called externally.
       */
      template<typename lookup_table_t>
      inline double posmaker( gsl_rng * r,
			      const double & beg,
			      const double & end,
			      lookup_table_t * lookup) const
      {
	double pos = gsl_ran_flat(r,beg,end);
	while( lookup->find(pos) != lookup->end() )
	  {
	    pos = gsl_ran_flat(r,beg,end);
	  }
	lookup->insert(pos);
	return pos;
      }

      /*!
	Return a KTfwd::popgenmut

	\param r A gsl_rng
	\param nmu Neutral mutation rate (per gamete, per generation)
	\param smu Selected mutation rate (per gamete, per generation)
	\param lookup Lookup table for mutations
       */
      template<typename lookup_table_t>
      inline result_type make_mut(gsl_rng * r,
				  const double & nmu,
				  const double & smu,
				  unsigned generation,
				  lookup_table_t * lookup) const
      {
	bool is_neutral = (gsl_rng_uniform(r) < nmu/(nmu+smu)) ? true : false;
	if( is_neutral )
	  {
	    size_t region = gsl_ran_discrete(r,nlookup.get());
	    double pos = posmaker(r,nbeg[region],nend[region],lookup);
	    return KTfwd::popgenmut(pos,0.,0.,generation,1);
	  }
	size_t region = gsl_ran_discrete(r,slookup.get());
	double pos = posmaker(r,sbeg[region],send[region],lookup);
	return KTfwd::popgenmut(pos,shmodels[region].s(r),shmodels[region].h(r),generation,1);
      }
    };

    struct discrete_rec_model
    /*!
      Class allowing the simulation of discrete variation
      in recombination rates along a region.
     */
    {
      using result_type = double;
      std::vector<double> beg,end;
      KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr lookup;
      /*!
	\param __beg Region beginnings
	\param __end Region ends
	\param __weight Region weights
       */
      discrete_rec_model( const std::vector<double> & __beg,
			  const std::vector<double> & __end,
			  const std::vector<double> & __weight ) : beg(__beg),end(__end)
      {
	lookup = KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(gsl_ran_discrete_preproc(__weight.size(),&__weight[0]));
      }
      /*!
	Returns a position from a region that is chosen based on region weights.
       */
      inline result_type operator()(gsl_rng * r) const
      {
	size_t region = gsl_ran_discrete(r,lookup.get());
	return gsl_ran_flat(r,beg[region],end[region]);
      }
    };
  }
}
#endif
