#ifndef FWDPP_SUGAR_POPBASE_HPP
#define FWDPP_SUGAR_POPBASE_HPP

/*!
  Base class for a population object

  Introduced in 0.4.8 to reduce code redundancy.
*/

#include <type_traits>
#include <fwdpp/type_traits.hpp>

namespace KTfwd
{
  namespace sugar
  {
      /*!
	\brief Abstraction of what is needed to simulate a population using fwdpp
	using an individual-based sampler from fwdpp

	All that is missing is the mutation_type and the container types.

	See @ref md_md_sugar for rationale, etc.

	\ingroup sugar
      */
      template<typename mutation_type,
	       typename mcont,
	       typename gcont,
	       typename mvector,
	       typename ftvector,
	       typename lookup_table_type,
	       typename tag_type>
      class popbase
      {
	static_assert(typename KTfwd::traits::is_gamete_t<typename gcont::value_type>::type(),
		      "gcont::value_type must be a gamete type" );
	static_assert(typename KTfwd::traits::is_mutation_t<typename mcont::value_type>::type(),
		      "mcont::value_type must be a mutation type" );
      public:
	uint_t N;

	//! Dispatch tag for other parts of sugar layer
	using popmodel_t = tag_type;
	//! Mutation type
	using mutation_t = mutation_type;
	//! Gamete type
	using gamete_t = typename gcont::value_type;
	//! Mutation vec type
	using mcont_t = mcont;
	//! Mutation count vector type
	using mcount_t = std::vector<uint_t>;
	//! Gamete vec type
	using gcont_t = gcont;
	//! Lookup table type for recording mutation positions, etc.
	using lookup_table_t = lookup_table_type;
	//! container type for fixations
	using mvector_t = mvector;
	//! container type for fixation times
	using ftvector_t = ftvector;
	mcont mutations;
	/*!
	  Used to keep track of mutation frequencies.

	  Should have memory reserved externally,
	  based on some good guess.
	*/
	mcount_t mcounts;
	gcont gametes;

	/*!
	  Vectors for holding copies of pointers to mutations during recombination.
	  The requirement to declare these was introduced in fwdpp 0.3.3.

	  In previous versions of the library, vectors like this had to be allocated
	  for every crossover event for every generation.  The result was an excessive
	  number of requests for memory allocation.

	  Now, we create the vector once per replicate.  Further, we will reserve memory
	  here, to minimize reallocs, etc., within fwdpp.

	  Internally, fwdpp's job is to make sure that this vector is appropriately
	  and efficiently cleared, but only when needed.

	  \note: if not using the sugar features, you can create these vectors
	  only once per simulation...
	*/
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

	//! Constructor
	popbase( const uint_t & popsize,
		 typename gamete_t::mutation_container::size_type reserve_size = 100) : N(popsize),
											//No muts in the population
											mutations(mcont_t()),
											mcounts(mcount_t()),
											//The population contains a single gamete in 2N copies
											gametes(gcont(1,gamete_t(2*popsize))),
											neutral(typename gamete_t::mutation_container()),
											selected(typename gamete_t::mutation_container()),
											mut_lookup(lookup_table_type()),
											fixations(mvector()),
											fixation_times(ftvector())
	{
	  //This is a good number for reserving,
	  //allowing for extra allocations when recycling is doing its thing
	  gametes.reserve(4*N);
	  //Reserve memory
	  neutral.reserve(reserve_size);
	  selected.reserve(reserve_size);
	}

	popbase( const popbase & ) = default;
	popbase( popbase && ) = default;
	
	bool operator==( const popbase & rhs ) const
	{
	  return this->mutations == rhs.mutations &&
	    this->mcounts == rhs.mcounts &&
	    this->gametes == rhs.gametes &&
	    this->fixations == rhs.fixations &&
	    this->fixation_times == rhs.fixation_times;
	}
      
	//! Empty all the containers
	void clear()
	{
	  mutations.clear();
	  mcounts.clear();
	  gametes.clear();
	  mut_lookup.clear();
	  fixations.clear();
	  fixation_times.clear();
	}
      };
    }
}

#endif
