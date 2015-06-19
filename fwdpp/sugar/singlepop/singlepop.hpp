#ifndef __FWDPP_SUGAR_SINGLEPOP_SINGLEPOP_HPP__
#define __FWDPP_SUGAR_SINGLEPOP_SINGLEPOP_HPP__

/*
  A structure representing a single Wright-Fisher population.
  The user initizializes it with a population size, N
*/

#include <type_traits>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/sugar/serialization.hpp>
#include <fwdpp/sugar/poptypes/tags.hpp>

namespace KTfwd {
  namespace sugar {
    /*!
      \brief Abstraction of what is needed to simulate a single population
      using an individual-based sampler from fwdpp
      
      All that is missing is the mutation_type and the container types.
      
      Does not allow copy construction/assignment!
      
      See @ref md_md_sugar for rationale, etc.

      \ingroup sugar
    */
    template<typename mutation_type,
	     typename mlist,
	     typename glist,
	     typename dipvector,
	     typename mvector,
	     typename ftvector,
	     typename lookup_table_type>
    class singlepop
    {
      static_assert( std::is_same< typename glist::value_type,
		     KTfwd::gamete_base< typename mlist::value_type, mlist > >::value,
		     "glist::value_type must be same as KTfwd::gamete_base< typename mlist::value_type, mlist >" );
    public:
      unsigned N;

      //! Dispatch tag for other parts of sugar layer
      using popmodel_t = sugar::SINGLEPOP_TAG;
      //! Mutation type
      using mutation_t = mutation_type;
      //! Gamete type
      using gamete_t = typename glist::value_type;
      //! Diploid vector type
      using dipvector_t = dipvector;
      //! Diploid type (std::pair<glist_t::iterator,glist_t::iterator>)
      using diploid_t = typename dipvector_t::value_type;
      //! Mutation list type
      using mlist_t = mlist;
      //! Gamete list type
      using glist_t = glist;
      //! Lookup table type for recording mutation positions, etc.
      using lookup_table_t = lookup_table_type;
      //! container type for fixations
      using mvector_t = mvector;
      //! container type for fixation times
      using ftvector_t = ftvector;

      mlist mutations;
      glist gametes;
      dipvector_t diploids;

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
      /*! \brief vector<unsigned> records times when mutation_ts 
	were added to mut_lookup
      */
      ftvector fixation_times;

      //! Constructor
      singlepop( const unsigned & popsize,
		 typename gamete_t::mutation_container::size_type reserve_size = 100) : N(popsize),
											//No muts in the population
											mutations(mlist()),
											//The population contains a single gamete in 2N copies
											gametes(glist(1,gamete_t(2*popsize))),
											//All N diploids contain the only gamete in the pop
											diploids(dipvector_t(popsize,diploid_t(gametes.begin(),gametes.begin()))),
											neutral(typename gamete_t::mutation_container()),
											selected(typename gamete_t::mutation_container()),
											mut_lookup(lookup_table_type()),
											fixations(mvector()),
											fixation_times(ftvector())
      {
	//Reserve memory
	neutral.reserve(reserve_size);
	selected.reserve(reserve_size);
      }

      //! Deleted
      singlepop( singlepop & ) = delete;
      //! Deleted
      singlepop( const singlepop & ) = delete;
      //! Move constructor
      singlepop( singlepop &&  ) = default;
      
      //! Deleted
      singlepop & operator=(singlepop &) = delete;
      //! Deleted
      singlepop & operator=(const singlepop &) = delete;

      //! Empty all the containers
      void clear() 
      {
	mutations.clear();
	gametes.clear();
	diploids.clear();
	mut_lookup.clear();
	fixations.clear();
	fixation_times.clear();
      }
    };

    /*!
      \brief Abstraction of what is needed to simulate a single population
      using an individual-based sampler from fwdpp.
      
      Allows copy/construction/assignment via deep copy.

      See @ref md_md_sugar for rationale, etc.

      \ingroup sugar
    */
    template<typename mutation_type,
	     typename mwriter,
	     typename mreader,
	     typename mlist,
	     typename glist,
	     typename dipvector,
	     typename mvector,
	     typename ftvector,
	     typename lookup_table_type,
	     typename dip_writer_t = KTfwd::diploidIOplaceholder,
	     typename dip_reader_t = KTfwd::diploidIOplaceholder>
    class singlepop_serialized
    {
      static_assert( std::is_same< typename glist::value_type,
		     KTfwd::gamete_base< typename mlist::value_type, mlist > >::value,
		     "glist::value_type must be same as KTfwd::gamete_base< typename mlist::value_type, mlist >" );
    public:
      unsigned N;

      //! Dispatch tag for other parts of sugar layer
      using popmodel_t = sugar::SINGLEPOP_TAG;
      
      //! Mutation type
      using mutation_t = mutation_type;
      //! Gamete type
      using gamete_t = typename glist::value_type;
      //! Vector of diploids
      using dipvector_t = dipvector;
      //! Diploid type (std::pair<glist_t::iterator,glist_t::iterator>)
      using diploid_t = typename dipvector_t::value_type;
      //! mutation_t serialization type
      using mwriter_t = mwriter;
      //! mutation_t deserialization type
      using mreader_t = mreader;
      //! Mutation list type
      using mlist_t = mlist;
      //! Gamete list type
      using glist_t = glist;
      //! Lookup table type for recording mutation positions, etc.
      using lookup_table_t = lookup_table_type;
      //! Serialization type for writing diploid genotypes
      using diploid_reader_t = dip_reader_t;
      //! Serialization type for reading diploid genotypes
      using diploid_writer_t = dip_writer_t;
      //! container type for fixations
      using mvector_t = mvector;
      //! container type for fixation times
      using ftvector_t = ftvector;

      //Data types -- the names should make the above typedefs a bit more clear
      mlist mutations;
      glist gametes;
      dipvector diploids;

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
      /*! \brief vector<unsigned> records times when mutation_ts 
	were added to mut_lookup
      */
      ftvector fixation_times;
      
      //! Constructor
      singlepop_serialized( const unsigned & popsize,
			    typename gamete_t::mutation_container::size_type reserve_size = 100) : N(popsize),
												   //No muts in the population
												   mutations(mlist()),
												   //The population contains a single gamete in 2N copies
												   gametes(glist(1,gamete_t(2*popsize))),
												   //All N diploids contain the only gamete in the pop
												   diploids(dipvector_t(popsize,diploid_t(gametes.begin(),gametes.begin()))),
												   neutral(typename gamete_t::mutation_container()),
												   selected(typename gamete_t::mutation_container()),
												   mut_lookup(lookup_table_type()),
												   fixations(mvector()),
												   fixation_times(ftvector())
      {
	neutral.reserve(reserve_size);
	selected.reserve(reserve_size);
      }
    
      //! Copy constructor
      singlepop_serialized(const singlepop_serialized & pop) : N(0u),
							       //No muts in the population
							       mutations(mlist()),
							       //The population contains a single gamete in 2N copies
							       gametes(glist(1,gamete_t(1))),
							       //All N diploids contain the only gamete in the pop
							       diploids(dipvector_t()),
							       neutral(pop.neutral),
							       selected(pop.selected),
							       mut_lookup(lookup_table_type()),
							       fixations(mvector()),
							       fixation_times(ftvector())
      {
	static_assert( std::is_same<mutation_t,typename mreader_t::result_type>::value,
		       "Mutation type must be same for class and mreader_t" );
	serialize s;
	s(pop,mwriter_t(),dip_writer_t());
	deserialize()(*this,s,mreader_t(),dip_reader_t());
      }
      
      //! Move constructor
      singlepop_serialized( singlepop_serialized && ) = default;

      //! Assignment operator
      singlepop_serialized & operator=(const singlepop_serialized & p)
      {
	static_assert( std::is_same<mutation_t,typename mreader_t::result_type>::value,
		       "Mutation type must be same for class and mreader_t" );
	serialize s;
	s(p,mwriter_t(),dip_writer_t());
	deserialize()(*this,s,mreader_t(),dip_reader_t());
	return *this;
      }

      //! Clear all containers
      void clear() 
      {
	mutations.clear();
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
