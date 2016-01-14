#ifndef _FITNESS_MODELS_HPP_
#define _FITNESS_MODELS_HPP_

#include <fwdpp/forward_types.hpp>
#include <fwdpp/fwd_functional.hpp>
#include <fwdpp/tags/diploid_tags.hpp>
#include <fwdpp/type_traits.hpp>
#include <cassert>
#include <type_traits>
#include <algorithm>
#include <functional>

/*!
  \defgroup fitness Policies for calculating fitnesses.
  This group contains the following data structures to help you implement custom fitness policies:
  1. KTfwd::site_dependent_fitness for implementing typical population-genetic-like models where trait values are a function of the properties of individual mutations
  2. KTfwd::haplotype_dependent_fitness for implementing models where trait values are the functions properties of haplotypes.

  Truthfully, latter is so trivial that a library user may never see any need for it.

  The library also defines two site-dependent fitness models:
  1. KTfwd::multiplicative_diploid
  2. KTfwd::additive_diploid

  These are arguably the "standard" models of the field as far as selection is concerned.

  Finally, KTfwd::no_selection is provided to force all diploid fitnesses to be equal to 1.
*/
namespace KTfwd
{
  /*! \brief Returns a fitness of 1
    \return A fitness of 1
    \param g1 A gamete
    \param g2 A gamete

    \note g1 and g2 must be part of the gamete_base hierarchy
    \ingroup fitness
  */
  struct no_selection
  {
    /*!
      \brief Method for standard diploid simulations of a single locus.
    */
    using result_type = double;
    template<typename gamete_type,
	     typename mcont_t>
    inline result_type operator()(const gamete_type &, const gamete_type &,
				  const mcont_t & ) const noexcept
    {
      static_assert( traits::is_gamete_t<gamete_type>::value,
                     "gamete_type::value_type must be a gamete type" );
      static_assert( traits::is_mutation_t<typename mcont_t::value_type>::value,
		     "mcont_t::value_type must be a mutation type" );
      return 1.;
    }
    //! \brief Naive implementation for non-standard cases
    template<typename T >
    inline result_type operator()(const T &) const noexcept
    {
      return 1.;
    }
  };

  /*! \brief Function object for fitness as a function of individual mutations in a diploid

    Function object for fitness as a function of mutations in a diploid.  Examples include the standard multiplicative and additive models of population genetics.  This routine idenfifies all homozygous/heterozygous mutations in a diploid and updates the diploid's fitness according to user-defined policies.  See the code for KTfwd::multiplicative_diploid and KTfwd::additive_diploid for specific examples.
    \param g1 An gamete
    \param g2 An gamete
    \param mutations Container of mutations
    \param fpol_hom Policy for updating fitness for the case of homozygosity for a mutant
    \param fpol_het Policy for updating fitness for the case of heterozygosity for a mutant
    \param starting_fitness The value to which the function will initialize the return value
    \return The fitness of a diploid with genotype g1 and g2
    \note The updating policies must take a non-const reference to a double as the first argument and
    an mcont_t::value_type as the second.  Further, they must not return anything. Any remaining arguments needed should be passed via a mechanism such as std::bind and a function object, or via a lambda expression.  See KTfwd::multiplicative_diploid for an example implementation.
    \ingroup fitness
  */
  struct site_dependent_fitness
  {
    using result_type = double;
    ///\example diploid_fixed_sh_ind.cc
    template<typename gamete_type,
	     typename mcont_t,
	     typename fitness_updating_policy_hom,
	     typename fitness_updating_policy_het>
    inline result_type operator()( const gamete_type & g1,
				   const gamete_type & g2,
				   const mcont_t & mutations,
				   const fitness_updating_policy_hom & fpol_hom,
				   const fitness_updating_policy_het & fpol_het,
				   const double & starting_fitness  = 1. ) const noexcept
    {
      static_assert( traits::is_gamete_t<gamete_type>::value,
                     "gamete_type::value_type must be a gamete type" );
      static_assert( traits::is_mutation_t<typename mcont_t::value_type>::value,
		     "mcont_t::value_type must be a mutation type" );
      static_assert( std::is_convertible<
		     fitness_updating_policy_hom,
		     std::function<void(double &,const typename mcont_t::value_type &)>
		     >::value,
		     "decltype(fpol_hom) must be convertible to std::function<void(double &,const typename mcont_t::value_type" );
      static_assert( std::is_convertible<
		     fitness_updating_policy_het,
		     std::function<void(double &,const typename mcont_t::value_type &)>
		     >::value,
		     "decltype(fpol_het) must be convertible to std::function<void(double &,const typename mcont_t::value_type" );
      result_type w = starting_fitness;
      if( g1.smutations.empty() && g2.smutations.empty() ) return w;
      else if (g1.smutations.empty()) {
	for( const auto & i : g2.smutations ) fpol_het( w,mutations[i] );
	return w;
      } else if (g2.smutations.empty()) {
	for( const auto & i : g1.smutations ) fpol_het( w,mutations[i] );
	return w;
      }
      typename gamete_type::mutation_container::size_type b2 = 0;
      for(const auto & b1 : g1.smutations)
	{
	  bool found = false;
	  for( ; !found && b2 < g2.smutations.size() && !( mutations[g2.smutations[b2]].pos > mutations[b1].pos ) ; ++b2 )
	    {
	      if (b1==g2.smutations[b2])
		{
		  assert(mutations[b1].pos == mutations[g2.smutations[b2]].pos);
		  fpol_hom(w,mutations[b1]);
		  found=true;
		}
	      else
		{
		  assert(mutations[g2.smutations[b2]].pos < mutations[b1].pos);
		  fpol_het(w,mutations[g2.smutations[b2]]);
		}
	    }
	  if(!found) fpol_het(w,mutations[b1]);
	}
      for( ; b2 < g2.smutations.size() ; ++b2 ) fpol_het(w,mutations[g2.smutations[b2]]);
      
      return w;
    }
    
    
    /*!
      \brief Overload for custom diploids.  This is what a programmer's functions will call.
      \note See @ref md_md_customdip
    */
    template< typename diploid2dispatch,
	      typename fitness_updating_policy_hom,
	      typename fitness_updating_policy_het>
    inline result_type operator()( const diploid2dispatch & dip,
				   const fitness_updating_policy_hom & fpol_hom,
				   const fitness_updating_policy_het & fpol_het,
				   const double & starting_fitness = 1. ) const noexcept
    {
      return this->operator()(dip.first,dip.second,fpol_hom,fpol_het,starting_fitness);
    }
  };

  /*! \brief Function object for fitness as a function of the two haplotypes in a diploid

    \param g1 A gamete
    \param g2 A gamete
    \param mutation Container of mutations
    \param hpol A policy whose first argument is an iterator to a gamete. Remaining arguments may be bound via std::bind or the equivalent.  The policy returns a double representing the effect of this haplotype on fitness
    \param dpol A policy whose first two arguments are doubles, each of which represents the effect of g1 and g2, respectively.  Remaining arguments may be bound via std::bind or the equivalent.  The policy returns a double representing the fitness of a diploid genotype g1/g2
    \return dpol( hpol(g1), hpol(g2) )
    \note This really is just a convenience function. Depending on the specifics of the model, this function may be totally unnecessary.
    \ingroup fitness
  */
  struct haplotype_dependent_fitness
  {
    using result_type = double;
    template< typename gamete_type,
	      typename mcont_t,
	      typename haplotype_policy,
	      typename diploid_policy >
    inline result_type operator()(const gamete_type & g1,
				  const gamete_type & g2,
				  const mcont_t & mutations,
				  const haplotype_policy & hpol,
				  const diploid_policy & dpol) const noexcept
    {
      static_assert( typename traits::is_gamete_t<gamete_type>::type(),
                     "gamete_type must be a gamete type" );
      static_assert( traits::is_mutation_t<typename mcont_t::value_type>::value,
		     "mcont_t::value_type must be a mutation type" );
      return dpol( hpol(g1), hpol(g2) );
    }
  };

  /*! \brief Multiplicative fitness across sites
    \param g1 A gamete
    \param g2 A gamete
    \param mutation Container of mutations
    \param scaling Fitnesses are 1, 1+h*s, 1+scaling*s, for AA,Aa,aa, resp.  This parameter lets you make sure your
    simulation is on the same scale as various formula in the literature
    \return Multiplicative fitness across sites.
    \ingroup fitness
  */
  struct multiplicative_diploid
  {
    using result_type = double;
    template< typename gamete_type, typename mcont_t>
    inline double operator()(const gamete_type & g1, const gamete_type & g2,
			     const mcont_t & mutations,
			     const double & scaling = 1.) const noexcept
    {
      using __mtype =  typename mcont_t::value_type;
      return std::max(0.,site_dependent_fitness()(g1,g2,mutations,
						  [&](double & fitness,const __mtype & mut) noexcept
						  {
						    fitness *= (1. + scaling*mut.s);
						  },
						  [](double & fitness,const __mtype & mut) noexcept
						  {
						    fitness *= (1. + mut.h*mut.s);
						  },
						  1.));
    }
    /*!
      \brief Overload for custom diploids.  This is what a programmer's functions will call.
      \note See @ref md_md_customdip
    */
    template< typename diploid2dispatch,
	      typename gcont_t,typename mcont_t>
    inline result_type operator()( const diploid2dispatch & dip,
				   const gcont_t & gametes,
				   const mcont_t & mutations,
				   const double & scaling = 1. ) const noexcept
    {
      return this->operator()(gametes[dip.first],gametes[dip.second],mutations,scaling);
    }
  };

  /*! \brief Additive fitness across sites
    \param g1 A gamete
    \param g2 A gamete
    \param mutations A container of mutations
    \param scaling Fitnesses are 1, 1+h*s, 1+scaling*s, for AA,Aa,aa, resp.  This parameter lets you make sure your
    simulation is on the same scale as various formula in the literature
    \return Additive fitness across sites.
    \note g1 and g2 must be part of the gamete_base hierarchy
    \ingroup fitness
  */
  struct additive_diploid
  {
    using result_type = double;
    template< typename gamete_type,typename mcont_t>
    inline result_type operator()(const gamete_type & g1, const gamete_type & g2,
				  const mcont_t & mutations,
				  const double & scaling = 1.) const noexcept
    {
      static_assert( typename traits::is_gamete_t<gamete_type>::type(),
                     "gamete_type must be a gamete type" );
      static_assert( traits::is_mutation_t<typename mcont_t::value_type>::value,
		     "mcont_t::value_type must be a mutation type" );
      using __mtype =  typename mcont_t::value_type;
      return std::max(0.,1. + site_dependent_fitness()(g1,g2,mutations,
						       [=](double & fitness,const __mtype & mut) noexcept
						       {
							 fitness += (scaling*mut.s);
						       },
						       [](double & fitness,const __mtype & mut) noexcept
						       {
							 fitness += (mut.h*mut.s);
						       },
						       0.)
		      );
    }

    /*!
      \brief Overload for custom diploids.  This is what a programmer's functions will call.
      \note See @ref md_md_customdip
    */
    template< typename diploid2dispatch,
	      typename gcont_t,
	      typename mcont_t>
    inline result_type operator()( const diploid2dispatch & dip,
				   const gcont_t & gametes,
				   const mcont_t & mutations,
				   const double & scaling = 1. ) const noexcept
    {
      return this->operator()(gametes[dip.first],gametes[dip.second],scaling);
    }
  };
}
#endif /* _FITNESS_MODELS_HPP_ */
