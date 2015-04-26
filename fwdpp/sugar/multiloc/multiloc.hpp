#ifndef __FWDPP_SUGAR_MULTILOC_MULTILOC_HPP__
#define __FWDPP_SUGAR_MULTILOC_MULTILOC_HPP__

#include <type_traits>
#include <vector>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/sugar/serialization.hpp>
#include <fwdpp/sugar/poptypes/tags.hpp>
#include <iostream>
namespace KTfwd {
  namespace sugar {
    /*!
      \brief Abstraction of what is needed to simulate a multilocus
      simulation using an individual-based sampler from fwdpp
      
      All that is missing is the mutation_type and the container types.
      
      Does not allow copy construction/assignment!
    */
    template<typename mutation_type,
	     typename mlist,
	     typename glist,
	     typename dipvector,
	     typename vglist,
	     typename mvector,
	     typename ftvector,
	     typename lookup_table_type>
    struct multiloc
    {
      //Dispatch tags for other parts of sugar layer
      using popmodel_t = sugar::MULTILOCPOP_TAG;
      
      //Typedefs for various container
      using mutation_t = mutation_type;
      using gamete_t = typename glist::value_type;
      using dipvector_t = dipvector;
      using diploid_t = typename dipvector_t::value_type;
      using mlist_t = mlist;
      using glist_t = glist;
      using vglist_t = vglist;
      using lookup_table_t = lookup_table_type;
      
      //public, non-const data
      unsigned N;
      mlist_t mutations;
      vglist_t gametes;
      dipvector_t diploids;
      lookup_table_type mut_lookup;
      mvector fixations;
      ftvector fixation_times;

      multiloc(const unsigned & __N, const unsigned & __nloci) : N(__N),
								 mutations(mlist_t()),
								 gametes(vglist_t(__nloci,glist_t(1,gamete_t(2*__N)))),
								 diploids(dipvector_t()),
								 mut_lookup(lookup_table_t()),
								 fixations(mvector()),
								 fixation_times(ftvector())
      {
	diploid_t idip;
	for( auto gitr = gametes.begin() ; gitr != gametes.end() ; ++gitr )
	  {
	    idip.emplace_back(std::make_pair(gitr->begin(),gitr->begin()));
	  }
	diploids = dipvector_t(N,idip);
      }

      multiloc(multiloc &&) = default;

      multiloc( multiloc & ) = delete;
      multiloc( const multiloc & ) = delete;
      multiloc & operator=(multiloc &) = delete;
      multiloc & operator=(const multiloc &) = delete;

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
      \brief Abstraction of what is needed to simulate a multilocus
      simulation using an individual-based sampler from fwdpp
      
      Allows copy/construction/assignment via deep copy
    */
    template<typename mutation_type,
	     typename mwriter,
	     typename mreader,
	     typename mlist,
	     typename glist,
	     typename dipvector,
	     typename vglist,
	     typename mvector,
	     typename ftvector,
	     typename lookup_table_type>
    struct multiloc_serialized
    {
      //Dispatch tags for other parts of sugar layer
      using popmodel_t = sugar::MULTILOCPOP_TAG;
      
      //Typedefs for various container
      using mutation_t = mutation_type;
      using gamete_t = typename glist::value_type;
      using dipvector_t = dipvector;
      using diploid_t = typename dipvector_t::value_type;
      using mlist_t = mlist;
      using glist_t = glist;
      using vglist_t = vglist;
      using lookup_table_t = lookup_table_type;
      using mwriter_t = mwriter;
      using mreader_t = mreader;
      
      //public, non-const data
      unsigned N;
      mlist_t mutations;
      vglist_t gametes;
      dipvector_t diploids;
      lookup_table_type mut_lookup;
      mvector fixations;
      ftvector fixation_times;

      multiloc_serialized(const unsigned & __N, const unsigned & __nloci) : N(__N),
									    mutations(mlist_t()),
									    gametes(vglist_t(__nloci,glist_t(1,gamete_t(2*__N)))),
									    diploids(dipvector_t()),
									    mut_lookup(lookup_table_t()),
									    fixations(mvector()),
									    fixation_times(ftvector())
      {
	diploid_t idip;
	for( auto gitr = gametes.begin() ; gitr != gametes.end() ; ++gitr )
	  {
	    idip.emplace_back(std::make_pair(gitr->begin(),gitr->begin()));
	  }
	diploids = dipvector_t(N,idip);
      }

      multiloc_serialized(multiloc_serialized &&) = default;

      multiloc_serialized( const multiloc_serialized & __m ) : N(0u), mutations(mlist_t()),
							       gametes(vglist_t()),
							       diploids(dipvector_t()),
							       mut_lookup(lookup_table_t()),
							       fixations(mvector()),
							       fixation_times(ftvector())
      {
	serialize s;
	s(__m,mwriter_t());
	deserialize()(*this,s,mreader_t());
      }
      
      multiloc_serialized & operator=(const multiloc_serialized & __m)
      {
	serialize s;
	s(__m,mwriter_t());
	deserialize()(*this,s,mreader_t());
	return *this;
      }

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
