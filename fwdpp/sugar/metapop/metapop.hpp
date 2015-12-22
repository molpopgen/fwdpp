#ifndef __FWDPP_SUGAR_METAPOP_METAPOP_HPP__
#define __FWDPP_SUGAR_METAPOP_METAPOP_HPP__

#include <type_traits>
#include <vector>
#include <numeric>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/type_traits.hpp>
#include <fwdpp/sugar/serialization.hpp>
#include <fwdpp/sugar/poptypes/tags.hpp>

namespace KTfwd {
  namespace sugar {
    /*!
      \brief Abstraction of what is needed to simulate a metapopulation
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
	     typename vdipvector,
	     typename mvector,
	     typename ftvector,
	     typename lookup_table_type>
    class metapop
    {
      static_assert(typename KTfwd::traits::is_gamete_t<typename glist::value_type>::type(),
		    "glist::value_type must be a gamete type" );
      static_assert(typename KTfwd::traits::is_mutation_t<typename mlist::value_type>::type(),
		    "mlist::value_type must be a mutation type" );

    private:
      void init_vectors()
      {
	uint_t metapopsize = std::accumulate(Ns.begin(),Ns.end(),0);
	gametes.emplace_back( gamete_t(2*metapopsize) );
	auto gam = gametes.begin();
	for(uint_t i = 0 ; i < Ns.size() ; ++i )
	  {
	    diploids.emplace_back( dipvector_t(Ns[i],typename dipvector::value_type(gam,gam)) );
	  }
      }
    public:
      //! Dispatch tags for other parts of sugar layer
      using popmodel_t = sugar::METAPOP_TAG;
      

      //! Mutation type
      using mutation_t = mutation_type;
      //! Gamete type
      using gamete_t = typename glist::value_type;
      //! Diploid vector type
      using dipvector_t = dipvector;
      //! Vector of dipvectors (vector of diploids for each deme)
      using vdipvector_t = vdipvector;
      //! Diploid type (std::pair<glist_t::iterator,glist_t::iterator>)
      using diploid_t = typename dipvector_t::value_type;
      //! Mutation list type
      using mlist_t = mlist;
      //! Gamete list type
      using glist_t = glist;
      //! Container of glist_t (container of gametes lists for each deme)
      //! Lookup table type for recording mutation positions, etc.
      using lookup_table_t = lookup_table_type;
      //! container type for fixations
      using mvector_t = mvector;
      //! container type for fixation times
      using ftvector_t = ftvector;

      //! Deme sizes
      std::vector<uint_t> Ns;
      mlist_t mutations;
      glist_t gametes;
      vdipvector_t diploids;

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
      
      //! Construct with a list of deme sizes
      metapop( std::initializer_list<uint_t> __Ns,
	       typename gamete_t::mutation_container::size_type reserve_size = 100) :
	Ns(__Ns),
	mutations(mlist_t()),
	gametes(glist_t()),
	diploids(vdipvector_t()),
	neutral(typename gamete_t::mutation_container()),
	selected(typename gamete_t::mutation_container()),
	mut_lookup(lookup_table_type()),
	fixations(mvector()),
	fixation_times(ftvector())
      {
	init_vectors();
	neutral.reserve(reserve_size);
	selected.reserve(reserve_size);
      }

      //! Construct with array of deme sizes
      metapop(const uint_t * __Ns, const size_t num_Ns,
	      typename gamete_t::mutation_container::size_type reserve_size = 100) :
	Ns(std::vector<uint_t>(__Ns,__Ns+num_Ns)),
	mutations(mlist_t()),
	gametes(glist_t()),
	diploids(vdipvector_t()),
	neutral(typename gamete_t::mutation_container()),
	selected(typename gamete_t::mutation_container()),
	mut_lookup(lookup_table_type()),
	fixations(mvector()),
	fixation_times(ftvector())
      {
	//Ns.assign(__Ns,__Ns+num_Ns);
	init_vectors();
	neutral.reserve(reserve_size);
	selected.reserve(reserve_size);
      }

      //! Move constructor
      metapop( metapop && __m ) = default;

      //! Deleted
      metapop( metapop & ) = delete;
      //! Deleted
      metapop( const metapop & ) = delete;
      //! Deleted
      metapop & operator=(metapop &) = delete;
      //! Deleted
      metapop & operator=(const metapop &) = delete;

      //! Empty all containers
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
      \brief Abstraction of what is needed to simulate a metapopulation
      using an individual-based sampler from fwdpp
      
      Provides for copy construction and assignment via deep copy.

      See @ref md_md_sugar for rationale, etc.

      \ingroup sugar
    */
    template<typename mutation_type,
	     typename mwriter,
	     typename mreader,
	     typename mlist,
	     typename glist,
	     typename dipvector,
	     typename vdipvector,
	     typename mvector,
	     typename ftvector,
	     typename lookup_table_type,
	     typename dip_writer_t = KTfwd::diploidIOplaceholder,
	     typename dip_reader_t = KTfwd::diploidIOplaceholder>
    class metapop_serialized
    {
      static_assert(typename KTfwd::traits::is_gamete_t<typename glist::value_type>::type(),
		    "glist::value_type must be a gamete type" );
      static_assert(typename KTfwd::traits::is_mutation_t<typename mlist::value_type>::type(),
		    "mlist::value_type must be a mutation type" );
    private:
      void init_vectors()
      {
	uint_t metapopsize = std::accumulate(Ns.begin(),Ns.end(),0);
	gametes.emplace_back( gamete_t(2*metapopsize) );
	auto gam = gametes.begin();
	for(uint_t i = 0 ; i < Ns.size() ; ++i )
	  {
	    diploids.emplace_back( dipvector_t(Ns[i],typename dipvector::value_type(gam,gam)) );
	  }
      }
    public:
      //! Dispatch tags for other parts of sugar layer
      using popmodel_t = sugar::METAPOP_TAG;
      
      //! Mutation type
      using mutation_t = mutation_type;
      //! mutation_t serialization type
      using mwriter_t = mwriter;
      //! mutation_t deserialization type
      using mreader_t = mreader;
      //! Gamete type
      using gamete_t = typename glist::value_type;
      //! Diploid vector type
      using dipvector_t = dipvector;
      //! Vector of dipvectors (vector of diploids for each deme)
      using vdipvector_t = vdipvector;
      //! Diploid type (std::pair<glist_t::iterator,glist_t::iterator>)
      using diploid_t = typename dipvector_t::value_type;
      //! Mutation list type
      using mlist_t = mlist;
      //! Gamete list type
      using glist_t = glist;
      //! Container of glist_t (container of gametes lists for each deme)
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

      //! Deme sizes
      std::vector<uint_t> Ns;
      mlist_t mutations;
      glist_t gametes;
      vdipvector_t diploids;

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
      
      //! Construct with a list of population sizes
      metapop_serialized( std::initializer_list<uint_t> __Ns,
			  typename gamete_t::mutation_container::size_type reserve_size = 100) : Ns(__Ns),
												 mutations(mlist_t()),
												 gametes(glist_t()),
												 diploids(vdipvector_t()),
												 neutral(typename gamete_t::mutation_container()),
												 selected(typename gamete_t::mutation_container()),
												 mut_lookup(lookup_table_type()),
												 fixations(mvector()),
												 fixation_times(ftvector())
      {
	init_vectors();
	neutral.reserve(reserve_size);
	selected.reserve(reserve_size);
      }

      metapop_serialized(const uint_t * __Ns, const size_t num_Ns,
			 typename gamete_t::mutation_container::size_type reserve_size = 100) :
	Ns(std::vector<uint_t>(__Ns,__Ns+num_Ns)),
	mutations(mlist_t()),
	gametes(glist_t()),
	diploids(vdipvector_t()),
	neutral(typename gamete_t::mutation_container()),
	selected(typename gamete_t::mutation_container()),
	mut_lookup(lookup_table_type()),
	fixations(mvector()),
	fixation_times(ftvector())
      {
	//	Ns.assign(__Ns,__Ns+num_Ns);
	init_vectors();
	neutral.reserve(reserve_size);
	selected.reserve(reserve_size);
      }

      //! Move constructor
      metapop_serialized( metapop_serialized && ) = default;

      //! Copy constructor
      metapop_serialized( const metapop_serialized & __m) : Ns(std::vector<uint_t>()),
							    mutations(mlist_t()),
							    gametes(glist_t()),
							    diploids(vdipvector_t()),
							    mut_lookup(lookup_table_type()),
							    fixations(mvector()),
							    fixation_times(ftvector())
      {
	serialize s;
	s(__m,mwriter_t(),dip_writer_t());
	deserialize()(*this,s,mreader_t(),dip_reader_t());
      }
      
      //! Assignment operator
      metapop_serialized & operator=(const metapop_serialized & __m)
      {
	serialize s;
	s(__m,mwriter_t(),dip_writer_t());
	deserialize()(*this,s,mreader_t(),dip_reader_t());
	return *this;
      }

      //! Empty all containers
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
