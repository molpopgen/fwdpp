#ifndef __KTFWD_MIGRATION_HPP__
#define __KTFWD_MIGRATION_HPP__

#include <fwdpp/forward_types.hpp>
#include <fwdpp/util.hpp>
#include <type_traits>
#include <algorithm>

namespace KTfwd
{

  /*!
    Migrate a randomly-chosen gamete from source to destination

    A gamete in destination is chosen at random and replaced
    with a gamete from source.

    \param r a GSL random number generator
    \param source A vector of gametes
    \param dest A vector of gametes
    \param twoN_1  2N for source
    \param twoN_1  2N for dest
    \param nmigrants The number of gametes to migrate from source to dest

    \example diploid_twopop_mig.cc
  */
  template<typename gamete_type,
	   typename vector_type_allocator,
	   template< typename,typename > class vector_type >
  void migrate_from_to( gsl_rng * r,
			vector_type<gamete_type,vector_type_allocator> * source, 	   
			vector_type<gamete_type,vector_type_allocator> * dest, 
			const unsigned & nmigrants)
  {
    typedef gamete_base< typename gamete_type::mutation_type, typename gamete_type::mutation_list_type > gamete_base_type;
    static_assert( std::is_base_of<gamete_base_type,gamete_type>::value ||
                   std::is_same<gamete_base_type,gamete_type >::value,
                   "gamete_type must be, or inherit from, KTfwd::gamete_base<mutation_type,mutation_list_type>" );
    typedef typename vector_type<gamete_type,vector_type_allocator>::iterator vgitr;
    for(unsigned i = 0 ; i < nmigrants ; ++i)
      {
	vgitr itr =  KTfwd::pgam( r, source ); 
	vgitr itr2 =  KTfwd::pgam( r, dest ); 
       	itr2->n--;
	gamete_type ng(*itr);
	ng.n=1;
	vgitr itr3 = update_if_exists_insert(ng,dest);  
	assert(std::count(dest->begin(),dest->end(),*itr3)==1);
      }
  }

  /*!
    Symmetric migration between two populations at rate m
    per generation.

    Here, m is defined such that Poissomn(twoN*m) gametes
    migrate from pop1 to pop2 and then from pop2 to pop1.

    A migration event is implemented via KTfwd::migrate_from_to,
    which defines a migration event as the replacement of a gamete in
    the destiniation population with a gamete from the source.  The frequency 
    of the source population gamete remains unchanged.

    \param r a GSL random number generator
    \param gametes1 Pointer to gametes from source population
    \param gametes2 Pointer to gametes from destination population
    \param twoN_1  2N for gametes1
    \param twoN_2  2N for gametes2
    \param m probability that a chromosome in population B is replaced by one from population A
   */
  template<typename gamete_type,
	   typename vector_type_allocator,
	   template< typename,typename > class vector_type >
  void migrate( gsl_rng * r,
		vector_type<gamete_type,vector_type_allocator> * gametes1, 	   
		vector_type<gamete_type,vector_type_allocator> * gametes2, 
		const unsigned & twoN_1, const unsigned & twoN_2,
		const double & m)
  {
    typedef gamete_base< typename gamete_type::mutation_type, typename gamete_type::mutation_list_type > gamete_base_type;
    static_assert( std::is_base_of<gamete_base_type,gamete_type>::value ||
                   std::is_same<gamete_base_type,gamete_type >::value,
                   "gamete_type must be, or inherit from, KTfwd::gamete_base<mutation_type,mutation_list_type>" );
    unsigned m12 = gsl_ran_poisson(r,double(twoN_1)*m);
    unsigned m21 = gsl_ran_poisson(r,double(twoN_2)*m);

    migrate_from_to(r, gametes1, gametes2, m12 );
    migrate_from_to(r, gametes2, gametes1, m21 );
  }
}


#endif
