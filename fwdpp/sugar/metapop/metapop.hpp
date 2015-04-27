#ifndef __FWDPP_SUGAR_METAPOP_METAPOP_HPP__
#define __FWDPP_SUGAR_METAPOP_METAPOP_HPP__

#include <type_traits>
#include <vector>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/sugar/serialization.hpp>
#include <fwdpp/sugar/poptypes/tags.hpp>

namespace KTfwd {
  namespace sugar {
    /*!
      \brief Abstraction of what is needed to simulate a metapopulation
      using an individual-based sampler from fwdpp
      
      All that is missing is the mutation_type and the container types.
      
      Does not allow copy construction/assignment!
    */
    template<typename mutation_type,
	     typename mlist,
	     typename glist,
	     typename dipvector,
	     typename vglist,
	     typename vdipvector,
	     typename mvector,
	     typename ftvector,
	     typename lookup_table_type>
    class metapop
    {
      static_assert( std::is_same< typename glist::value_type,
		     KTfwd::gamete_base< typename mlist::value_type, mlist > >::value,
		     "glist::value_type must be same as KTfwd::gamete_base< typename mlist::value_type, mlist >" );
    private:
      void init_vectors()
      {
	for( unsigned i = 0 ; i < Ns.size() ; ++i )
	  {
	    gametes.emplace_back( glist_t(1,gamete_t(2*Ns[i])) );
	    diploids.emplace_back(dipvector_t(Ns[i],std::make_pair( gametes[i].begin(),gametes[i].begin())));
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
      //! Diploid type (std::pair<glist_t::iterator,glist_t::iterator)
      using diploid_t = typename dipvector_t::value_type;
      //! Mutation list type
      using mlist_t = mlist;
      //! Gamete list type
      using glist_t = glist;
      //! Container of glist_t (container of gametes lists for each deme)
      using vglist_t = vglist;
      //! Lookup table type for recording mutation positions, etc.
      using lookup_table_t = lookup_table_type;

      //! Deme sizes
      std::vector<unsigned> Ns;
      mlist_t mutations;
      vglist_t gametes;
      vdipvector_t diploids;
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
      
      //! Construct with a list of deme sizes
      metapop( std::initializer_list<unsigned> __Ns ) : Ns(__Ns),
							mutations(mlist_t()),
							gametes(vglist_t()),
							diploids(vdipvector_t()),
							mut_lookup(lookup_table_type()),
							fixations(mvector()),
							fixation_times(ftvector())
      {
	init_vectors();
      }

      //! Construct with array of deme sizes
      metapop(const unsigned * __Ns, const size_t num_Ns) : Ns(std::vector<unsigned>()),
							    mutations(mlist_t()),
							    gametes(vglist_t()),
							    diploids(vdipvector_t()),
							    mut_lookup(lookup_table_type()),
							    fixations(mvector()),
							    fixation_times(ftvector())
      {
	Ns.assign(__Ns,__Ns+num_Ns);
	init_vectors();
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
    */
    template<typename mutation_type,
	     typename mwriter,
	     typename mreader,
	     typename mlist,
	     typename glist,
	     typename dipvector,
	     typename vglist,
	     typename vdipvector,
	     typename mvector,
	     typename ftvector,
	     typename lookup_table_type>
    class metapop_serialized
    {
      static_assert( std::is_same< typename glist::value_type,
		     KTfwd::gamete_base< typename mlist::value_type, mlist > >::value,
		     "glist::value_type must be same as KTfwd::gamete_base< typename mlist::value_type, mlist >" );
    private:
      void init_vectors()
      {
	for( unsigned i = 0 ; i < Ns.size() ; ++i )
	  {
	    gametes.emplace_back( glist_t(1,gamete_t(2*Ns[i])) );
	    diploids.emplace_back(dipvector_t(Ns[i],std::make_pair( gametes[i].begin(),gametes[i].begin())));
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
      //! Diploid type (std::pair<glist_t::iterator,glist_t::iterator)
      using diploid_t = typename dipvector_t::value_type;
      //! Mutation list type
      using mlist_t = mlist;
      //! Gamete list type
      using glist_t = glist;
      //! Container of glist_t (container of gametes lists for each deme)
      using vglist_t = vglist;
      //! Lookup table type for recording mutation positions, etc.
      using lookup_table_t = lookup_table_type;

      //! Deme sizes
      std::vector<unsigned> Ns;
      mlist_t mutations;
      vglist_t gametes;
      vdipvector_t diploids;
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
      
      //! Construct with a list of population sizes
      metapop_serialized( std::initializer_list<unsigned> __Ns ) : Ns(__Ns),
								   mutations(mlist_t()),
								   gametes(vglist_t()),
								   diploids(vdipvector_t()),
								   mut_lookup(lookup_table_type()),
								   fixations(mvector()),
								   fixation_times(ftvector())
      {
	init_vectors();
      }

      metapop_serialized(const unsigned * __Ns, const size_t num_Ns) : Ns(std::vector<unsigned>()),
								       mutations(mlist_t()),
								       gametes(vglist_t()),
								       diploids(vdipvector_t()),
								       mut_lookup(lookup_table_type()),
								       fixations(mvector()),
								       fixation_times(ftvector())
      {
	Ns.assign(__Ns,__Ns+num_Ns);
	init_vectors();
      }

      //! Move constructor
      metapop_serialized( metapop_serialized && ) = default;

      //! Copy constructor
      metapop_serialized( const metapop_serialized & __m) : Ns(std::vector<unsigned>()),
							    mutations(mlist_t()),
							    gametes(vglist_t()),
							    diploids(vdipvector_t()),
							    mut_lookup(lookup_table_type()),
							    fixations(mvector()),
							    fixation_times(ftvector())
      {
	serialize s;
	s(__m,mwriter_t());
	deserialize()(*this,s,mreader_t());
      }
      
      //! Assignment operator
      metapop_serialized & operator=(const metapop_serialized & __m)
      {
	serialize s;
	s(__m,mwriter_t());
	deserialize()(*this,s,mreader_t());
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
