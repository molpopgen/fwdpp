#ifndef __FWDPP_EXTENSIONS_REGIONS_HPP__
#define __FWDPP_EXTENSIONS_REGIONS_HPP__

#include <fwdpp/sugar/infsites.hpp>
#include <fwdpp/internal/gsl_discrete.hpp>
#include <fwdpp/extensions/callbacks.hpp>
#include <iostream>
namespace KTfwd
{
  namespace extensions
  {
    struct discrete_mut_model
    {
      using result_type = KTfwd::popgenmut;
      std::vector<double> nbeg,nend,sbeg,send;
      std::vector< shmodel > shmodels;
      KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr nlookup,slookup;
      discrete_mut_model( const std::vector<double> &__nbeg,
			  const std::vector<double> &__nend,
			  const std::vector<double> & nweights, //the weights
			  const std::vector<double> &__sbeg,
			  const std::vector<double> &__send,
			  const std::vector<double> & sweights, //the weights
			  const std::vector<shmodel> &__shmodels) : nbeg(__nbeg),nend(__nend),
								   sbeg(__sbeg),send(__send),
								   shmodels(__shmodels)
      {
	std::cerr << "fwdpp here1\n";
	std::cerr << nbeg.size() << ' ' << nend.size() << ' ' << nweights.size() << '\n';
	if(nweights.size())
	  nlookup = KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(gsl_ran_discrete_preproc(nweights.size(),&nweights[0]));
	std::cerr << "fwdpp here2\n";
	if(sweights.size())
	  slookup = KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(gsl_ran_discrete_preproc(sweights.size(),&sweights[0]));
      }

      //Helper fnx
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

      template<typename lookup_table_t>
      inline result_type make_mut(gsl_rng * r,
				  const double & nmu,
				  const double & smu,
				  const unsigned & generation,
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
    {
      using result_type = double;
      std::vector<double> beg,end;
      KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr lookup;
      discrete_rec_model( const std::vector<double> & __beg,
			  const std::vector<double> & __end,
			  const std::vector<double> & __weight ) : beg(__beg),end(__end)
      {
	lookup = KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(gsl_ran_discrete_preproc(__weight.size(),&__weight[0]));
      }
      inline result_type operator()(gsl_rng * r) const
      {
	size_t region = gsl_ran_discrete(r,lookup.get());
	return gsl_ran_flat(r,beg[region],end[region]);
      }
    };
  }
}
#endif
