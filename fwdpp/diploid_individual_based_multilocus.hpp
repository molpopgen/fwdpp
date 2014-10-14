#ifndef __FWDPP_IND_BASED_MLOC_HPP__
#define __FWDPP_IND_BASED_MLOC_HPP__

#include <utility>

namespace KTfwd
{
//single deme, N changing
template< typename gamete_type,
	  typename glist_vector_type_allocator,
	  typename gamete_list_type_allocator,
	  typename mutation_list_type_allocator,
	  typename diploid_vector_type_allocator,
	  typename locus_vector_type_allocator,
	  typename diploid_fitness_function,
	  typename mutation_removal_policy,
	  typename mutation_model_container,
	  typename recombination_policy_container,
	  typename mutation_insertion_policy,
	  typename gamete_insertion_policy,
	  template<typename,typename> class gamete_list_type,
	  template<typename,typename> class glist_vector_type,
	  template<typename,typename> class mutation_list_type,
	  template<typename,typename> class diploid_vector_type,
	  template<typename,typename> class locus_vector_type>
double
sample_diploid(gsl_rng * r,
	       glist_vector_type< gamete_list_type<gamete_type,
						   gamete_list_type_allocator> ,
	       glist_vector_type_allocator > * gametes,
	       diploid_vector_type<locus_vector_type<std::pair<typename gamete_list_type< gamete_type,gamete_list_type_allocator >::iterator,
							       typename gamete_list_type< gamete_type,gamete_list_type_allocator >::iterator>,
						     locus_vector_type_allocator>,
	       diploid_vector_type_allocator> * diploids,
	       mutation_list_type<typename gamete_type::mutation_type,mutation_list_type_allocator > * mutations, 
	       const unsigned & N_curr, 
	       const unsigned & N_next, 
	       const double * mu,
	       const mutation_model_container & mmodel,
	       const recombination_policy_container & rec_policies,
	       const double * r_between_loci,
	       const mutation_insertion_policy & mpolicy,
	       const gamete_insertion_policy & gpolicy_mut,
	       const diploid_fitness_function & ff,
	       const mutation_removal_policy & mp,
	       const double & f);

  //single deme, constant N
  template< typename gamete_type,
	    typename glist_vector_type_allocator,
	    typename gamete_list_type_allocator,
	    typename mutation_list_type_allocator,
	    typename diploid_vector_type_allocator,
	    typename locus_vector_type_allocator,
	    typename diploid_fitness_function,
	    typename mutation_removal_policy,
	    typename mutation_model_container,
	    typename recombination_policy_container,
	    typename mutation_insertion_policy,
	    typename gamete_insertion_policy,
	    template<typename,typename> class gamete_list_type,
	    template<typename,typename> class glist_vector_type,
	    template<typename,typename> class mutation_list_type,
	    template<typename,typename> class diploid_vector_type,
	    template<typename,typename> class locus_vector_type>
  double
  sample_diploid(gsl_rng * r,
		 glist_vector_type< gamete_list_type<gamete_type,
		 gamete_list_type_allocator> ,
		 glist_vector_type_allocator > * gametes,
		 diploid_vector_type<locus_vector_type<std::pair<typename gamete_list_type< gamete_type,gamete_list_type_allocator >::iterator,
		 typename gamete_list_type< gamete_type,gamete_list_type_allocator >::iterator>,
						       locus_vector_type_allocator>,
				     diploid_vector_type_allocator> * diploids,
		 mutation_list_type<typename gamete_type::mutation_type,mutation_list_type_allocator > * mutations, 
		 const unsigned & N,
		 const double * mu,
		 const mutation_model_container & mmodel,
		 const recombination_policy_container & rec_policies,
		 const double * r_between_loci,
		 const mutation_insertion_policy & mpolicy,
		 const gamete_insertion_policy & gpolicy_mut,
		 const diploid_fitness_function & ff,
		 const mutation_removal_policy & mp,
		 const double & f);
}
#endif
