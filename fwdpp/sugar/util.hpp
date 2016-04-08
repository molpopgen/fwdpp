/*!
  \file fwdpp/sugar/util.hpp

  \brief Header for miscellaneous functions 
*/

#ifndef FWDPP_SUGAR_UTIL_HPP
#define FWDPP_SUGAR_UTIL_HPP

#include <exception>
#include <type_traits>
#include <algorithm>
#include <fwdpp/internal/recycling.hpp>
#include <fwdpp/sugar/poptypes/tags.hpp>

namespace KTfwd
{
  namespace fwdpp_internal
  {
    template<typename poptype,
	     typename dipvector_t>
    /*!
      This function def'n looks funny b/c we use 
      diploids to refer both to the diploids in a single-deme sim
      or the diploids in a particular deme of a meta-pop sim.

      \note Do not call directly.  Use KTfwd::add_mutation instead.
    */
    void add_mutation_details( poptype & p,
			       dipvector_t & diploids,
			       typename poptype::mcont_t::value_type & new_mutation,
			       const std::vector<std::size_t> & indlist,
			       const std::vector<short> & clist)
    {
      //find an extinct mutation, if one exists
      auto extinct_mut = std::find(std::begin(p.mcounts),std::end(p.mcounts),0);
      auto gam_recycling_bin = fwdpp_internal::make_gamete_queue(p.gametes);

      //Add new mutation to mutation container and track its location
      std::size_t mindex=p.mcounts.size();  //this will be the correct value if we use the else block below
      if( extinct_mut != std::end(p.mcounts) ) //then we can recycle
	{
	  auto dist=std::distance(std::begin(p.mcounts),extinct_mut);
	  p.mutations[dist] = std::move(new_mutation); //move new mutation into place
	  mindex=std::size_t(dist); //update our value
	}
      else
	{
	  //cannot recycle, so add it to end
	  p.mutations.emplace_back(std::move(new_mutation));
	  p.mcounts.push_back(0); //Add a place for this variant
	}

      //Function object for calls to upper bound
      auto inserter = [&p](const double & __value,const std::size_t __mut) noexcept {
	assert(__mut<p.mutations.size());
	return __value < p.mutations[__mut].pos;
      };
      
      //update the diploids
      bool neutral = p.mutations[mindex].neutral;
      auto pos = p.mutations[mindex].pos;
      for( std::size_t i = 0 ; i < indlist.size() ; ++i )
	{
	  if( clist[i] == 0 || clist[i] == 2 )
	    {
	      auto n = p.gametes[diploids[indlist[i]].first].mutations;
	      auto s = p.gametes[diploids[indlist[i]].first].smutations;
	      if(neutral)
		{
		  n.insert( std::upper_bound(n.begin(),
					     n.end(),pos,
					     inserter),
			    mindex );
		}
	      else
		{
		  s.insert( std::upper_bound(s.begin(),
					     s.end(),pos,
					     inserter),
			    mindex );
		}
	      diploids[indlist[i]].first = fwdpp_internal::recycle_gamete(p.gametes,
									gam_recycling_bin,
									n,s);
	      p.mcounts[mindex]++;
	    }
	  if(clist[i]>0)
	    {
	      auto n = p.gametes[diploids[indlist[i]].second].mutations;
	      auto s = p.gametes[diploids[indlist[i]].second].smutations;
	      if(p.mutations[mindex].neutral)
		{
		  n.insert( std::upper_bound(n.begin(),
					     n.end(),pos,
					     inserter),
			    mindex );
		}
	      else
		{
		  s.insert( std::upper_bound(s.begin(),
					     s.end(),pos,
					     inserter),
			    mindex );
		}
	      diploids[indlist[i]].second = fwdpp_internal::recycle_gamete(p.gametes,
									   gam_recycling_bin,
									   n,s);
	      p.mcounts[mindex]++;
	    }
	}
    }      
  }
  
  template<typename poptype,
	   class... Args>
  void add_mutation(poptype & p,
		    const std::vector<std::size_t> & indlist,
		    const std::vector<short> & clist,
		    Args&&... args)
  /*!
    \brief Add a mutation into a population at a given frequency.

    \param p A single deme object.
    \param indlist A list of indexes of diploids into which to add the new mutations.
    \param clist A list of gametes.  See below.
    \param args Values required to cosnstruct a new mutation.  See below.

    Some notes:

    clist.size() must equal indlist.size()

    Values in \a clist must be 0, 1, or 2. These values mean to add the mutation to the first,
    second, or both gametes, resepectively, of each diploid in \a indlist.

    See the unit test file unit/test_sugar_add_mutation.cc for example of use.
  */
  {
    static_assert( std::is_same<typename poptype::popmodel_t,KTfwd::sugar::SINGLEPOP_TAG>::value,
		   "poptype must be a single-deme object type" );

    //Before we go deep into creating objects, let's do some checks
    for(const auto & i : indlist)
      {
	if(i>=p.diploids.size()) throw std::out_of_range("indlist contains elements > p.diploids.size()");
      }
    for( const auto & c : clist )
      {
	if(c<0||c>2) throw std::runtime_error("clist contains elements < 0 and/or > 2");
      }
    if(indlist.size()!=clist.size()) throw std::runtime_error("indlist and clist must be same length");
    
    //create a new mutation
    typename poptype::mcont_t::value_type new_mutant(std::forward<Args>(args)...);
    fwdpp_internal::add_mutation_details(p,p.diploids,new_mutant,indlist,clist);
  }
}

#endif
