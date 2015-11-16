#ifndef __FWDPP_SUGAR_MUTATION_INFSITES_HPP__
#define __FWDPP_SUGAR_MUTATION_INFSITES_HPP__

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fwdpp/sugar/popgenmut.hpp>

namespace KTfwd
{
  /*!
    \brief Generic function object implementing 
    the infinitely-many sites mutation model for 
    "standard" population-genetic scenarios
    \ingroup sugar
   */
  struct infsites
  {
    /*!
      \param r A gsl_rng *
      \param lookup A lookup table of mutation positions, see @ref md_md_policies for 
      \param generation Generation when this mutation is happening
      \param neutral_mutation_rate Either the rate at which neutral variants arise (per gamete per generation), or something directly proportional to it
      \param selected_mutation_rate Either the rate at which non-neutral variants arise (per gamete per generation), or something directly proportional to it
      \param posmaker A policy that returns the position of the new mutation
      \param smaker A policy generating the selection coefficient/effect size associated with non-neutral variants
      \param hmaker A policy generating the dominance associated with non-neutral variants

      \note A mutation will be "selected" with probability selected_mutation_rate/(selected_mutation_rate + neutral_mutation_rate)
     */
    template<typename lookup_table_t,
	     typename position_t,
	     typename sdist_t,
	     typename hdist_t>
    inline popgenmut 
    operator()(gsl_rng * r, lookup_table_t * lookup,
	       const uint_t & generation,
	       const double & neutral_mutation_rate,
	       const double & selected_mutation_rate,
	       const position_t & posmaker,
	       const sdist_t & smaker,
	       const hdist_t & hmaker) const
    {
      //Establish position of new mutation
      double pos = posmaker();
      while(lookup->find(pos) != lookup->end())
	{
	  pos = posmaker();
	}
      lookup->insert(pos);
      //Is mutation selected or not?
      if( gsl_rng_uniform(r) < selected_mutation_rate/(neutral_mutation_rate + selected_mutation_rate) )
	{
	  return popgenmut(pos,smaker(),hmaker(),generation,1);
	}
      //return a neutral mutation
      return popgenmut(pos,0.,0.,generation,1);
    }

    /*!
      \brief Overload for different position distributions for neutral and non-neutral variants

      \param r A gsl_rng *
      \param lookup A lookup table of mutation positions, see @ref md_md_policies for 
      \param generation Generation when this mutation is happening
      \param neutral_mutation_rate Either the rate at which neutral variants arise (per gamete per generation), or something directly proportional to it
      \param selected_mutation_rate Either the rate at which non-neutral variants arise (per gamete per generation), or something directly proportional to it
      \param nposmaker A policy that returns the position of new neutral mutation
      \param sposmaker A policy that returns the position of new non-neutral mutation
      \param smaker A policy generating the selection coefficient/effect size associated with non-neutral variants
      \param hmaker A policy generating the dominance associated with non-neutral variants

      \note A mutation will be "selected" with probability selected_mutation_rate/(selected_mutation_rate + neutral_mutation_rate)
     */
    template<typename lookup_table_t,
	     typename nposition_t,
	     typename sposition_t,
	     typename sdist_t,
	     typename hdist_t>
    inline popgenmut
    operator()(gsl_rng * r, lookup_table_t * lookup,
	       const uint_t & generation,
	       const double & neutral_mutation_rate,
	       const double & selected_mutation_rate,
	       const nposition_t & nposmaker,
	       const sposition_t & sposmaker,
	       const sdist_t & smaker,
	       const hdist_t & hmaker) const
    {
      bool selected = gsl_rng_uniform(r) < selected_mutation_rate/(neutral_mutation_rate + selected_mutation_rate);
      if( selected )
	{
	  double pos = sposmaker();
	  while(lookup->find(pos) != lookup->end())
	    {
	      pos = sposmaker();
	    }
	  lookup->insert(pos);
	  return popgenmut(pos,smaker(),hmaker(),generation,1);
	}
      //Establish position of new mutation
      double pos = nposmaker();
      while(lookup->find(pos) != lookup->end())
	{
	  pos = nposmaker();
	}
      lookup->insert(pos);
      //return a neutral mutation
      return popgenmut(pos,0.,0.,generation,1);
    }

    /*!
      \param r A gsl_rng *
      \param lookup A lookup table of mutation positions, see @ref md_md_policies for 
      \param generation Generation when this mutation is happening
      \param neutral_mutation_rate Either the rate at which neutral variants arise (per gamete per generation), or something directly proportional to it
      \param selected_mutation_rate Either the rate at which non-neutral variants arise (per gamete per generation), or something directly proportional to it
      \param posmaker A policy that returns the position of the new mutation
      \param smaker A policy generating the selection coefficient/effect size associated with non-neutral variants
      \param hmaker A policy generating the dominance associated with non-neutral variants

      \note A mutation will be "selected" with probability selected_mutation_rate/(selected_mutation_rate + neutral_mutation_rate)
     */
    template<typename lookup_table_t,
	     typename position_t,
	     typename sdist_t,
	     typename hdist_t>
    inline popgenmut
    operator()(gsl_rng * r, lookup_table_t * lookup,
	       const uint_t * generation,
	       const double & neutral_mutation_rate,
	       const double & selected_mutation_rate,
	       const position_t & posmaker,
	       const sdist_t & smaker,
	       const hdist_t & hmaker) const
    {
      //Establish position of new mutation
      double pos = posmaker();
      while(lookup->find(pos) != lookup->end())
	{
	  pos = posmaker();
	}
      lookup->insert(pos);
      //Is mutation selected or not?
      if( gsl_rng_uniform(r) < selected_mutation_rate/(neutral_mutation_rate + selected_mutation_rate) )
	{
	  return popgenmut(pos,smaker(),hmaker(),*generation,1);
	}
      //return a neutral mutation
      return popgenmut(pos,0.,0.,*generation,1);
    }

    /*!
      \param r A gsl_rng *
      \param lookup A lookup table of mutation positions, see @ref md_md_policies for 
      \param neutral_mutation_rate Either the rate at which neutral variants arise (per gamete per generation), or something directly proportional to it
      \param selected_mutation_rate Either the rate at which non-neutral variants arise (per gamete per generation), or something directly proportional to it
      \param posmaker A policy that returns the position of the new mutation
      \param smaker A policy generating the selection coefficient/effect size associated with non-neutral variants
      \param hmaker A policy generating the dominance associated with non-neutral variants

      \note A mutation will be "selected" with probability selected_mutation_rate/(selected_mutation_rate + neutral_mutation_rate)
     */
    template<typename lookup_table_t,
	     typename position_t,
	     typename sdist_t,
	     typename hdist_t>
    inline mutation
    operator()(gsl_rng * r, lookup_table_t * lookup,
	       const double & neutral_mutation_rate,
	       const double & selected_mutation_rate,
	       const position_t & posmaker,
	       const sdist_t & smaker,
	       const hdist_t & hmaker) const
    {
      //Establish position of new mutation
      double pos = posmaker();
      while(lookup->find(pos) != lookup->end())
	{
	  pos = posmaker();
	}
      lookup->insert(pos);
      //Is mutation selected or not?
      if( gsl_rng_uniform(r) < selected_mutation_rate/(neutral_mutation_rate + selected_mutation_rate) )
	{
	  return mutation(pos,smaker(),1,hmaker());
	}
      //return a neutral mutation
      return mutation(pos,0.,1,0.);
    }

    /*!
      \brief Overload for different position distributions for neutral and non-neutral variants

      \param r A gsl_rng *
      \param lookup A lookup table of mutation positions, see @ref md_md_policies for 
      \param neutral_mutation_rate Either the rate at which neutral variants arise (per gamete per generation), or something directly proportional to it
      \param selected_mutation_rate Either the rate at which non-neutral variants arise (per gamete per generation), or something directly proportional to it
      \param nposmaker A policy that returns the position of new neutral mutation
      \param sposmaker A policy that returns the position of new non-neutral mutation
      \param smaker A policy generating the selection coefficient/effect size associated with non-neutral variants
      \param hmaker A policy generating the dominance associated with non-neutral variants

      \note A mutation will be "selected" with probability selected_mutation_rate/(selected_mutation_rate + neutral_mutation_rate)
    */
    template<typename lookup_table_t,
	     typename nposition_t,
	     typename sposition_t,
	     typename sdist_t,
	     typename hdist_t>
    inline mutation
    operator()(gsl_rng * r, lookup_table_t * lookup,
	       const double & neutral_mutation_rate,
	       const double & selected_mutation_rate,
	       const nposition_t & nposmaker,
	       const sposition_t & sposmaker,
	       const sdist_t & smaker,
	       const hdist_t & hmaker) const
    {
      bool selected = gsl_rng_uniform(r) < selected_mutation_rate/(neutral_mutation_rate + selected_mutation_rate);
      if ( selected )
	{
	  double pos = sposmaker();
	  while(lookup->find(pos) != lookup->end())
	    {
	      pos = sposmaker();
	    }
	  lookup->insert(pos);
	  return mutation(pos,smaker(),1,hmaker());
	}
      double pos = nposmaker();
      while(lookup->find(pos) != lookup->end())
	{
	  pos = nposmaker();
	}
      lookup->insert(pos);
      //return a neutral mutation
      return mutation(pos,0.,1,0.);
    }
  };
}

#endif
