#ifndef __FWDPP_SUGAR_SINGLEPOP_HPP__
#define __FWDPP_SUGAR_SINGLEPOP_HPP__

/*
  A structure representing a single Wright-Fisher population.
  The user initizializes it with a population size, N
*/

#include <type_traits>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/sugar/serialization.hpp>

namespace KTfwd {
  namespace sugar {
    /*!
      Abstraction of what is needed to simulate a single population
      using an individual-based sampler from fwdpp

      All that is missing is the mutation_type and the container types.

      Does not allow copy construction/assignment!
    */
    template<typename mutation_type,
	     typename mlist,
	     typename glist,
	     typename mvector,
	     typename ftvector,
	     typename lookup_table_type>
    class singlepop
    {
    public:
      unsigned N;

      //Typedefs for various container
      using mtype = mutation_type;
      using gtype = KTfwd::gamete_base< typename mlist::value_type, mlist >;
      using dipvector = std::vector< std::pair<typename glist::iterator,typename glist::iterator> >;
      using mlist_t = mlist;
      using glist_t = glist;
      using lookup_table_t = lookup_table_type;
      
      //Data types -- the names should make the above typedefs a bit more clear
      mlist mutations;
      glist gametes;
      dipvector diploids;
      lookup_table_type mut_lookup;
      mvector fixations;
      ftvector fixation_times;

      //Constructors
      singlepop( const unsigned & popsize ) : N(popsize),
					      mutations(mlist()),                //No muts in the population
					      gametes(glist(1,gtype(2*popsize))), //The population contains a single gamete in 2N copies
					      diploids(dipvector(popsize,std::make_pair(gametes.begin(),gametes.begin()))), //All N diploids contain the only gamete in the pop
					      mut_lookup(lookup_table_type()),
					      fixations(mvector()),
					      fixation_times(ftvector())
      {
      }

      //Do NOT allow copy construction
      singlepop( singlepop & ) = delete;
      singlepop( const singlepop & ) = delete;
      //Allow move construction
      singlepop( singlepop && p ) : N (std::move(p.N)),
				    mutations(std::move(p.mutations)),
				    gametes(std::move(p.gametes)),
				    diploids(std::move(p.diploids)),
				    fixations(std::move(p.fixations)),
				    fixation_times(std::move(p.fixation_times))
      {
	//Fill the mutation lookup!
	std::for_each( mutations.begin(), mutations.end(),
		       [this]( const mtype & __m ) { mut_lookup.insert(__m.pos); } );
      }
      //Do not allow assignment from a reference
      singlepop & operator=(singlepop &) = delete;
      singlepop & operator=(const singlepop &) = delete;

      //Member functions
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
      Abstraction of what is needed to simulate a single population
      using an individual-based sampler from fwdpp
      
      All that is missing is the mutation_type and the container types.
      
      Allows copy/construction/assignment via deep copy
    */
    template<typename mutation_type,
	     typename mwriter,
	     typename mreader,
	     typename mlist,
	     typename glist,
	     typename mvector,
	     typename ftvector,
	     typename lookup_table_type>
    class singlepop_serialized
    {
    public:
      unsigned N;

      //Typedefs for various container
      using mtype = mutation_type;
      using gtype = KTfwd::gamete_base< typename mlist::value_type, mlist >;
      using dipvector = std::vector< std::pair<typename glist::iterator,typename glist::iterator> >;
      using mwriter_t = mwriter;
      using mreader_t = mreader;
      using mlist_t = mlist;
      using glist_t = glist;
      using lookup_table_t = lookup_table_type;
      
      //Data types -- the names should make the above typedefs a bit more clear
      mlist mutations;
      glist gametes;
      dipvector diploids;
      lookup_table_type mut_lookup;
      mvector fixations;
      ftvector fixation_times;
      
      //Constructors
      singlepop_serialized( const unsigned & popsize ) : N(popsize),
							 mutations(mlist()),                //No muts in the population
							 gametes(glist(1,gtype(2*popsize))), //The population contains a single gamete in 2N copies
							 diploids(dipvector(popsize,std::make_pair(gametes.begin(),gametes.begin()))), //All N diploids contain the only gamete in the pop
							 mut_lookup(lookup_table_type()),
							 fixations(mvector()),
							 fixation_times(ftvector())
      {
      }

      singlepop_serialized( const singlepop_serialized & pop) : N(0u),
								mutations(mlist()),                //No muts in the population
								gametes(glist(1,gtype(1))), //The population contains a single gamete in 2N copies
								diploids(dipvector()), //All N diploids contain the only gamete in the pop
								mut_lookup(lookup_table_type()),
								fixations(mvector()),
								fixation_times(ftvector())
      {
	static_assert( std::is_same<mtype,typename mreader_t::result_type>::value,
		       "Mutation type must be same for class and mreader_t" );
	serialize s;
	s(pop,mwriter_t());
	deserialize()(*this,s,mreader_t());
      }
      //Allow move construction
      singlepop_serialized( singlepop_serialized && p ) : N (std::move(p.N)),
				    mutations(std::move(p.mutations)),
				    gametes(std::move(p.gametes)),
				    diploids(std::move(p.diploids)),
				    fixations(std::move(p.fixations)),
				    fixation_times(std::move(p.fixation_times))
      {
	//Fill the mutation lookup!
	std::for_each( mutations.begin(), mutations.end(),
		       [this]( const mtype & __m ) { mut_lookup.insert(__m.pos); } );
      }

      singlepop_serialized & operator=(const singlepop_serialized & p)
      {
	static_assert( std::is_same<mtype,typename mreader_t::result_type>::value,
		       "Mutation type must be same for class and mreader_t" );
	serialize s;
	s(p,mwriter_t());
	deserialize()(*this,s,mreader_t());
	return *this;
      }

      //Member functions
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
