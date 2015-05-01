#ifndef __FWDPP_SUGAR_MUTATION_INFSITES_HPP__
#define __FWDPP_SUGAR_MUTATION_INFSITES_HPP__

#include <type_traits>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fwdpp/sugar/popgenmut.hpp>

namespace KTfwd
{
  struct infsites
  {
    template<typename mlist_t,
	     typename lookup_table_t,
	     typename position_t,
	     typename sdist_t,
	     typename hdist_t>
    inline typename std::enable_if<std::is_same<typename mlist_t::value_type,popgenmut>::value,typename mlist_t::value_type>::type
    operator()(gsl_rng * r, mlist_t * mutations, lookup_table_t * lookup,
	       const unsigned & generation,
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
      if( gsl_rng_uniform(r) <= selected_mutation_rate/(neutral_mutation_rate + selected_mutation_rate) )
	{
	  return typename mlist_t::value_type(pos,smaker(),hmaker(),generation,1);
	}
      //return a neutral mutation
      return typename mlist_t::value_type(pos,0.,0.,generation,1);
    }

    /*!
      Overload for different position generators for neutral and non-neutral variants
    */
    template<typename mlist_t,
	     typename lookup_table_t,
	     typename nposition_t,
	     typename sposition_t,
	     typename sdist_t,
	     typename hdist_t>
    inline typename std::enable_if<std::is_same<typename mlist_t::value_type,popgenmut>::value,typename mlist_t::value_type>::type
    operator()(gsl_rng * r, mlist_t * mutations, lookup_table_t * lookup,
	       const unsigned & generation,
	       const double & neutral_mutation_rate,
	       const double & selected_mutation_rate,
	       const nposition_t & nposmaker,
	       const sposition_t & sposmaker,
	       const sdist_t & smaker,
	       const hdist_t & hmaker) const
    {
      bool selected = gsl_rng_uniform(r) <= selected_mutation_rate/(neutral_mutation_rate + selected_mutation_rate);
      if( selected )
	{
	  double pos = sposmaker();
	  while(lookup->find(pos) != lookup->end())
	    {
	      pos = sposmaker();
	    }
	  lookup->insert(pos);
	  return typename mlist_t::value_type(pos,smaker(),hmaker(),generation,1);
	}
      //Establish position of new mutation
      double pos = nposmaker();
      while(lookup->find(pos) != lookup->end())
	{
	  pos = nposmaker();
	}
      lookup->insert(pos);
      //return a neutral mutation
      return typename mlist_t::value_type(pos,0.,0.,generation,1);
    }

    template<typename mlist_t,
	     typename lookup_table_t,
	     typename position_t,
	     typename sdist_t,
	     typename hdist_t>
    inline typename std::enable_if<std::is_same<typename mlist_t::value_type,popgenmut>::value,typename mlist_t::value_type>::type
    operator()(gsl_rng * r, mlist_t * mutations, lookup_table_t * lookup,
	       const unsigned * generation,
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
      if( gsl_rng_uniform(r) <= selected_mutation_rate/(neutral_mutation_rate + selected_mutation_rate) )
	{
	  return typename mlist_t::value_type(pos,smaker(),hmaker(),*generation,1);
	}
      //return a neutral mutation
      return typename mlist_t::value_type(pos,0.,0.,*generation,1);
    }
    
    template<typename mlist_t,
	     typename lookup_table_t,
	     typename position_t,
	     typename sdist_t,
	     typename hdist_t>
    inline typename std::enable_if<std::is_same<typename mlist_t::value_type,mutation>::value,typename mlist_t::value_type>::type
    operator()(gsl_rng * r, mlist_t * mutations, lookup_table_t * lookup,
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
      if( gsl_rng_uniform(r) <= selected_mutation_rate/(neutral_mutation_rate + selected_mutation_rate) )
	{
	  return typename mlist_t::value_type(pos,smaker(),1,hmaker());
	}
      //return a neutral mutation
      return typename mlist_t::value_type(pos,0.,1,0.);
    }

    /*!
      Overload for different position distributions for neutral and non-neutral variants
    */
    template<typename mlist_t,
	     typename lookup_table_t,
	     typename nposition_t,
	     typename sposition_t,
	     typename sdist_t,
	     typename hdist_t>
    inline typename std::enable_if<std::is_same<typename mlist_t::value_type,mutation>::value,typename mlist_t::value_type>::type
    operator()(gsl_rng * r, mlist_t * mutations, lookup_table_t * lookup,
	       const double & neutral_mutation_rate,
	       const double & selected_mutation_rate,
	       const nposition_t & nposmaker,
	       const sposition_t & sposmaker,
	       const sdist_t & smaker,
	       const hdist_t & hmaker) const
    {
      bool selected = gsl_rng_uniform(r) <= selected_mutation_rate/(neutral_mutation_rate + selected_mutation_rate);
      if ( selected )
	{
	  double pos = sposmaker();
	  while(lookup->find(pos) != lookup->end())
	    {
	      pos = sposmaker();
	    }
	  lookup->insert(pos);
	  return typename mlist_t::value_type(pos,smaker(),1,hmaker());
	}
      double pos = nposmaker();
      while(lookup->find(pos) != lookup->end())
	{
	  pos = nposmaker();
	}
      lookup->insert(pos);
      //return a neutral mutation
      return typename mlist_t::value_type(pos,0.,1,0.);
    }
  };
}

#endif
