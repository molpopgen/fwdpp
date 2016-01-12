#ifndef __FWDPP_SUGAR_MULTILOC_MULTILOC_HPP__
#define __FWDPP_SUGAR_MULTILOC_MULTILOC_HPP__

#include <type_traits>
#include <vector>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/sugar/poptypes/tags.hpp>

namespace KTfwd {
  namespace sugar {
    /*!
      \brief Abstraction of what is needed to simulate a multilocus
      simulation using an individual-based sampler from fwdpp.
      
      All that is missing is the mutation_type and the container types.
      
      Does not allow copy construction/assignment!

      See @ref md_md_sugar for rationale, etc.

      \ingroup sugar
    */
    template<typename mutation_type,
	     typename mcont,
	     typename gcont,
	     typename dipvector,
	     typename mvector,
	     typename ftvector,
	     typename lookup_table_type>
    struct multiloc
    {
      //! Dispatch tags for other parts of sugar layer
      using popmodel_t = sugar::MULTILOCPOP_TAG;
      
      //! Mutation type
      using mutation_t = mutation_type;
      //! Gamete type
      using gamete_t = typename gcont::value_type;
      //! Vector of vector of diploids
      //using vdipvector_t = vdipvector;
      //! Vector of diploids
      using dipvector_t = dipvector;
      //! Diploid type 
      using diploid_t = typename dipvector_t::value_type;
      //! Linked cont of mutations
      using mcont_t = mcont;
      //! Mutation count vector type
      using mcount_t = std::vector<uint_t>;
      //! Linked cont of gametes
      using gcont_t = gcont;
      //! Container of gcont_t (gametes for each locus)
      //using vgcont_t = vgcont;
      //! Lookup table type for recording mutation positions, etc.
      using lookup_table_t = lookup_table_type;
      //! container type for fixations
      using mvector_t = mvector;
      //! container type for fixation times
      using ftvector_t = ftvector;

      //! Population size
      uint_t N;
      mcont_t mutations;
      /*!
	Used to keep track of mutation frequencies.

	Should have memory reserved externally,
	based on some good guess.
      */
      mcount_t mcounts;
      gcont_t gametes;
      dipvector_t diploids;
      //! Vectors for recombination intermediates
      typename gamete_t::mutation_container neutral,selected;
      /*!
	\brief Can be used to track positions of segregating mutations.
	\note Must have interface like std::map or std::unordered_set
      */
      lookup_table_type mut_lookup;
      //! Vector of mutation_t to track fixations
      mvector fixations;
      /*! \brief vector<uint_t> records times when mutation_ts 
	were added to mut_lookup
      */
      ftvector fixation_times;

      //! Construct with population size and number of loci
      multiloc(const uint_t & __N, const uint_t & __nloci,
	       typename gamete_t::mutation_container::size_type reserve_size = 100) : N(__N),
										      mutations(mcont_t()),
										      mcounts(mcount_t()),
										      gametes(gcont_t(1,gamete_t(2*__N))),
										      diploids(dipvector_t()),
										      mut_lookup(lookup_table_t()),
										      fixations(mvector()),
										      fixation_times(ftvector())
      {
	diploid_t d(__nloci,typename diploid_t::value_type(0,0));
	diploids = dipvector_t(N,d);
	neutral.reserve(reserve_size);
	selected.reserve(reserve_size);
      }

      bool operator==( const metapop & rhs ) const
      {
	return this->mutations == rhs.mutations &&
	  this->mcounts == rhs.mcounts &&
	  this->gametes == rhs.gametes &&
	  this->diploids == rhs.diploids &&
	  this->fixations == rhs.fixations &&
	  this->fixation_times == rhs.fixation_times;
      }
      
      //! Empty all containers
      void clear() 
      {
	mutations.clear();
	mcounts.clear();
	gametes.clear();
	diploids.clear();
	mut_lookup.clear();
	fixations.clear();
	fixation_times.clear();
      }
    };

  }
}

#endif
