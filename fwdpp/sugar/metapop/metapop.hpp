#ifndef __FWDPP_SUGAR_METAPOP_METAPOP_HPP__
#define __FWDPP_SUGAR_METAPOP_METAPOP_HPP__

#include <type_traits>
#include <vector>
#include <numeric>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/type_traits.hpp>
#include <fwdpp/sugar/singlepop/singlepop.hpp>
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
	     typename mcont,
	     typename gcont,
	     typename dipvector,
	     typename vdipvector,
	     typename mvector,
	     typename ftvector,
	     typename lookup_table_type>
    class metapop
    {
      static_assert(typename KTfwd::traits::is_gamete_t<typename gcont::value_type>::type(),
		    "gcont::value_type must be a gamete type" );
      static_assert(typename KTfwd::traits::is_mutation_t<typename mcont::value_type>::type(),
		    "mcont::value_type must be a mutation type" );
    private:
      void init_vectors()
      {
	uint_t metapopsize = std::accumulate(Ns.begin(),Ns.end(),0u);
	gametes.emplace_back( gamete_t(2*metapopsize) );
	for(uint_t i = 0 ; i < Ns.size() ; ++i )
	  {
	    diploids.emplace_back( dipvector_t(Ns[i],typename dipvector::value_type(0,0)) );
	  }
      }
    public:
      //! Dispatch tags for other parts of sugar layer
      using popmodel_t = sugar::METAPOP_TAG;
      

      //! Mutation type
      using mutation_t = mutation_type;
      //! Gamete type
      using gamete_t = typename gcont::value_type;
      //! Diploid vector type
      using dipvector_t = dipvector;
      //! Vector of dipvectors (vector of diploids for each deme)
      using vdipvector_t = vdipvector;
      //! Diploid type (std::pair<gcont_t::iterator,gcont_t::iterator>)
      using diploid_t = typename dipvector_t::value_type;
      //! Mutation cont type
      using mcont_t = mcont;
      //! Mutation count vector type
      using mcount_t = std::vector<uint_t>;
      //! Gamete cont type
      using gcont_t = gcont;
      //! Lookup table type for recording mutation positions, etc.
      using lookup_table_t = lookup_table_type;
      //! container type for fixations
      using mvector_t = mvector;
      //! container type for fixation times
      using ftvector_t = ftvector;
      //! Fitness function signature compatible with this type
      using fitness_t = KTfwd::traits::fitness_fxn_t<dipvector_t,gcont_t,mcont_t>;
      //! Metapops can be constructed from singlepops of this type
      using compat_singlepop_t = sugar::singlepop<mutation_t,
						  mcont_t,
						  gcont_t,
						  dipvector_t,
						  mvector_t,
						  ftvector_t,
						  lookup_table_t>;
      //! Deme sizes
      std::vector<uint_t> Ns;
      mcont_t mutations;
      /*!
	Used to keep track of mutation frequencies.

	Should have memory reserved externally,
	based on some good guess.
      */
      mcount_t mcounts;
      gcont_t gametes;
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
      
      //! Construct with a cont of deme sizes
      metapop( std::initializer_list<uint_t> __Ns,
	       typename gamete_t::mutation_container::size_type reserve_size = 100) :
	Ns(__Ns),
	mutations(mcont_t()),
	mcounts(mcount_t()),
	gametes(gcont_t()),
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
	mutations(mcont_t()),
	mcounts(mcount_t()),
	gametes(gcont_t()),
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
						  
      //! Copy construct from a singlepop based on the same basic types
      metapop( const compat_singlepop_t & spop ) : Ns({spop.N}),
						   mutations(spop.mutations),
						   mcounts(spop.mcounts),
						   gametes(spop.gametes),
						   diploids(vdipvector_t(1,spop.diploids)),
						   neutral(spop.neutral),
						   selected(spop.selected),
						   mut_lookup(spop.mut_lookup),
						   fixations(spop.fixations),
						   fixation_times(spop.fixation_times)
      {
      }

      //! Move construct from a singlepop based on the same basic types
      metapop( compat_singlepop_t && spop ) : Ns({spop.N}),
					      mutations(std::move(spop.mutations)),
					      mcounts(std::move(spop.mcounts)),
					      gametes(std::move(spop.gametes)),
					      diploids(vdipvector_t(1,std::move(spop.diploids))),
					      neutral(std::move(spop.neutral)),
					      selected(std::move(spop.selected)),
					      mut_lookup(std::move(spop.mut_lookup)),
					      fixations(std::move(spop.fixations)),
					      fixation_times(std::move(spop.fixation_times))
      {
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
