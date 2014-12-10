#ifndef _FITNESS_MODELS_HPP_
#define _FITNESS_MODELS_HPP_

#include <fwdpp/forward_types.hpp>
#include <fwdpp/fwd_functional.hpp>
#include <cassert>
#include <type_traits>
#include <algorithm>
#include <functional>

namespace KTfwd
{
  /*! \brief Returns a fitness of 1
    \return A fitness of 1
    \param g1 A gamete
    \param g2 A gamete

    \note g1 and g2 must be part of the gamete_base hierarchy
  */
  struct no_selection
  {
    typedef double result_type;
    template<typename iterator_type >
    inline result_type operator()(const iterator_type & g1, const iterator_type &g2) const
    {
      static_assert( std::is_base_of<mutation_base,
		     typename iterator_type::value_type::mutation_type>::value,
                     "iterator_type::value_type::mutation_type must be derived from KTfwd::mutation_base" );
      return 1.;
    }
  };

  /*! \brief Function object for fitness as a function of individual mutations in a diploid

    Function object for fitness as a function of mutations in a diploid.  Examples include the standard multiplicative and additive models of population genetics.  This routine idenfifies all homozygous/heterozygous mutations in a diploid and updates the diploid's fitness according to user-defined policies.  See the code for KTfwd::multiplicative_diploid and KTfwd::additive_diploid for specific examples.
    \param g1 An iterator to a gamete
    \param g2 An iterator to a gamete
    \param fpol_hom Policy for updating fitness for the case of homozygosity for a mutant
    \param fpol_het Policy for updating fitness for the case of heterozygosity for a mutant
    \param starting_fitness The value to which the function will initialize the return value
    \return The fitness of a diploid with genotype g1 and g2
    \note The updating policies must take a non-const reference to a double as the first argument and
    an iterator to a gamete as the second.  Any remaining arguments needed should be passed via a mechanism such as std::bind and a function object, or via a lambda expression.  See KTfwd::multiplicative_diploid for an example implementation.
    \example diploid_fixed_sh_ind.cc
   */
  struct site_dependent_fitness
  {
    typedef double result_type;
    template<typename iterator_type,
	     typename fitness_updating_policy_hom,
	     typename fitness_updating_policy_het>
    inline result_type operator()( const iterator_type & g1,
				   const iterator_type & g2,
				   const fitness_updating_policy_hom & fpol_hom,
				   const fitness_updating_policy_het & fpol_het,
				   const double starting_fitness = 1. ) const
    {
      static_assert( std::is_base_of<mutation_base,
                                     typename iterator_type::value_type::mutation_type>::value,
                     "iterator_type::value_type::mutation_type must be derived from KTfwd::mutation_base" );
      result_type fitness=starting_fitness;
      if( g1->smutations.empty() && g2->smutations.empty() ) return fitness;
      if( g1->smutations.empty() && g2->smutations.empty() ) 
	{
	  std::for_each( g1->smutations.begin(),
			 g1->smutations.end(),
			 std::bind(fpol_het,std::ref(fitness),std::placeholders::_1) );
	  return fitness;
	}
      if( g1->smutations.empty() && !g2->smutations.empty() ) 
	{
	  std::for_each( g2->smutations.begin(),
			 g2->smutations.end(),
			 std::bind(fpol_het,std::ref(fitness),std::placeholders::_1) );
	  return fitness;
	}
      typename iterator_type::value_type::mutation_list_type_iterator ib1,ib2;
      typename iterator_type::value_type::mutation_container::const_iterator b1=g1->smutations.cbegin(),
	e1=g1->smutations.cend(),
	b2=g2->smutations.cbegin(),
	e2=g2->smutations.cend();
      //This is a fast way to calculate fitnesses,
      //as it just compares addresses in memory, and 
      //does little in the way of dereferencing and storing
      bool found = false;
      for( ; b1 < e1 ; ++b1 )
	{
	  found = false;
	  ib1 = *b1;
	  for( ; !found && b2 < e2 && !((*b2)->pos > (ib1)->pos) ; ++b2 )
	    {
	      ib2 = *b2;
	      if ( ib2 == ib1 ) //homozygote
		{
		  assert(ib1->s == ib2->s); //just making sure
		  assert(ib1->pos == ib2->pos);
		  fpol_hom(fitness,ib1);
		  found=true;
		}
	      else 
		//b2 points to a unique mutation that comes before b1
		{
		  assert(ib2->pos != ib1->pos);
		  assert(ib2->pos < ib1->pos);
		  fpol_het(fitness,ib2);
		}
	    }
	  if(!found) //het
	    {
	      fpol_het(fitness,ib1);
	    }
	}
      std::for_each( b2,e2,
		     std::bind(fpol_het,std::ref(fitness),std::placeholders::_1) );
      return fitness;
    }
  };

  template<typename iterator_type,
	   typename fitness_updating_policy_hom,
	   typename fitness_updating_policy_het>
    inline double site_dependent_fitness2( iterator_type __first1, iterator_type __last1,
					   iterator_type __first2, iterator_type __last2,
					   const fitness_updating_policy_hom &fpol_hom,
					   const fitness_updating_policy_het &fpol_het,
					   const double & starting_fitness = 1.)
  {
    double fitness = starting_fitness;
    // auto __x = [&](const typename iterator_type::value_type & mut_i)
    //   {
    // 	bool found = false;
    // 	while(! found && __first2 < __last2 )
    // 	  {
    // 	    if(mut_i == *__first2)
    // 	      {
    // 		found = true;
    // 		fpol_hom(fitness,mut_i);
    // 	      }
    // 	    else if( (*__first2)->pos < mut_i->pos )
    // 	      {
    // 		fpol_het(fitness,*__first2);
    // 	      }
    // 	    else
    // 	      {
    // 		break;
    // 	      }
    // 	    ++__first2;
    // 	  }
    // 	if(!found)
    // 	  {
    // 	    fpol_het(fitness,mut_i);
    // 	  }
    //   };
    // std::for_each( __first1, __last1,__x );
    
    while( __first1 < __last1 )
      {
    	bool found = false;
    	const double cpos = (*__first1)->pos;
    	while( !found && __first2 < __last2 ) 
    	  {
    	    if (*__first1 == *__first2) //homozygote
    	      {
    		found = true;
    		fpol_hom(fitness,*__first1);
    		//++__first2;
    		//break;
    	      }
    	    else if( (*__first2)->pos < cpos ) //then first2 points to a het
    	      {
    		fpol_het(fitness,*__first2);
    	      }
    	    else //first 1 is a het
    	      {
    		//fpol_het(fitness,*__first1);
    		//++__first2;
    		break;
    	      }
    	    ++__first2;
    	  }
    	if(!found) fpol_het(fitness,*__first1);
    	++__first1;
      }
    //std::cerr << std::distance(__first2,__last2) << '\n';
    // while(__first2<__last2) 
    //   {
    // 	fpol_het(fitness,*__first2);
    // 	++__first2;
    //   }
     std::for_each( __first2,__last2, 
		   std::bind(fpol_het,std::ref(fitness),std::placeholders::_1) );
     return fitness;
  }

  /*! \brief Function object for fitness as a function of the two haplotypes in a diploid

    \param g1 An iterator to a gamete
    \param g2 An iterator to a gamete
    \param hpol A policy whose first argument is an iterator to a gamete. Remaining arguments may be bound via std::bind or the equivalent.  The policy returns a double representing the effect of this haplotype on fitness
    \param dpol A policy whose first two arguments are doubles, each of which represents the effect of g1 and g2, respectively.  Remaining arguments may be bound via std::bind or the equivalent.  The policy returns a double representing the fitness of a diploid genotype g1/g2
    \return dpol( hpol(g1), hpol(g2) )
    \note This really is just a convenience function. Depending on the specifics of the model, this function may be totally unnecessary.
  */
  struct haplotype_dependent_fitness
  {
    typedef double result_type;
    template< typename iterator_type,
	      typename haplotype_policy,
	      typename diploid_policy >
    inline result_type operator()(const iterator_type & g1,
				  const iterator_type & g2,
				  const haplotype_policy & hpol,
				  const diploid_policy & dpol) const
    {
      static_assert( std::is_base_of<mutation_base,
                                     typename iterator_type::value_type::mutation_type>::value,
                     "iterator_type::value_type::mutation_type must be derived from KTfwd::mutation_base" );
      return dpol( hpol(g1), hpol(g2) );
    }
  };

  /*! \brief Multiplicative fitness across sites
    \param g1 An iterator pointing to a gamete
    \param g2 An iterator pointing to a gamete
    \param scaling Fitnesses are 1, 1+h*s, 1+scaling*s, for AA,Aa,aa, resp.  This parameter lets you make sure your
    simulation is on the same scale as various formula in the literature
    \return Multiplicative fitness across sites.
    1.);
  */
  struct multiplicative_diploid
  {
    typedef double result_type;
    template< typename iterator_type>
    inline double operator()(const iterator_type & g1, const iterator_type & g2, 
			     const double scaling = 1.) const
    {
      using __mtype =  typename iterator_type::value_type::mutation_list_type_iterator;
      return site_dependent_fitness()(g1,g2,
				      [&](double & fitness,const __mtype & mut)
				      {
					fitness *= (1. + scaling*mut->s);
				      },
				      [](double & fitness,const __mtype & mut)
				      {
					fitness *= (1. + mut->h*mut->s);
				      },
				      1.);
    }
  };

  /*! \brief Additive fitness across sites
    \param g1 An iterator pointing to a gamete
    \param g2 An iterator pointing to a gamete
    \param scaling Fitnesses are 1, 1+h*s, 1+scaling*s, for AA,Aa,aa, resp.  This parameter lets you make sure your
    simulation is on the same scale as various formula in the literature
    \return Additive fitness across sites.
    \note g1 and g2 must be part of the gamete_base hierarchy
  */
  struct additive_diploid
  {
    typedef double result_type;
    template< typename iterator_type>
    inline double operator()(const iterator_type & g1, const iterator_type & g2, 
			     const double scaling = 1.) const
    {
      using __mtype =  typename iterator_type::value_type::mutation_list_type_iterator;
      return 1. + site_dependent_fitness()(g1,g2,
					   [=](double & fitness,const __mtype & mut)
					   {
					     fitness += (scaling*mut->s);
					   },
					   [](double & fitness,const __mtype & mut)
					   {
					     fitness += (mut->h*mut->s);
					   },
					   0.);
    }
  };
}
#endif /* _FITNESS_MODELS_HPP_ */
