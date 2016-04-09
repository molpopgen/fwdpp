/*!
  \file fwdpp/sugar/util.hpp

  \brief Header for miscellaneous functions 
*/

#ifndef FWDPP_SUGAR_UTIL_HPP
#define FWDPP_SUGAR_UTIL_HPP

#include <exception>
#include <type_traits>
#include <algorithm>
#include <cassert>
#include <fwdpp/internal/recycling.hpp>
#include <fwdpp/sugar/poptypes/tags.hpp>

/*
  TODO:

  1. check if mutation to insert already exists.  Will this require operator== 
  to be defined...?  Done, and "yes"...

  2. re-work so that a vector of pointers to diploid gamete keys is collected and grouped by value,
  so that the minimum number of new gametes are created.  Ideally, the unit tests would change such that 
  only 1 new gamete is inserted for those examples.
*/

namespace KTfwd
{
  namespace sugar
  {
    template<typename mcont_t,
	     typename mcounts_t>
    std::size_t get_mut_index(mcont_t & mutations,
			      mcounts_t & mcounts,
			      typename mcont_t::value_type & new_mutation)
    /*!
      \brief Helper function for implementing KTfwd::add_mutation

      This function puts the new mutation into the mutations container and 
      updates other objects as needed.

      \return A key telling where \a new_mutation exists within \a mutations.

      \ingroup sugar
    */
    {
      //If the mutation already exists, let's not create it...
      auto mitr = std::find(std::begin(mutations),std::end(mutations),new_mutation);
      if(mitr!=mutations.end()) return std::size_t(std::distance(std::begin(mutations),mitr));

      /*
	If we are here, then new_mutation does not currently exist, so we have to add it.

	We do this using fwdpp's recycling method.  But, we don't need to create recycling bin,
	as we're just looking to add one mutation.  Instead, we can simply search for an extinct 
	mutation.
      */
      //find an extinct mutation, if one exists
      auto extinct_mut = std::find(std::begin(mcounts),std::end(mcounts),0);
      std::size_t mindex=mcounts.size();  //this will be the correct value if we use the else block below
      if( extinct_mut != std::end(mcounts) ) //then we can recycle
	{
	  auto dist=std::distance(std::begin(mcounts),extinct_mut);
	  mutations[dist] = std::move(new_mutation); //move new mutation into place
	  mindex=std::size_t(dist); //update our value
	}
      else
	{
	  //cannot recycle, so add it to end
	  mutations.emplace_back(std::move(new_mutation));
	  mcounts.push_back(0); //Add a place for this variant
	}
      return mindex;
    }
    
    template<typename poptype,
	     typename dipvector_t>
    /*!
      \brief Overload for single-deme sims and metapop sims.

      This function def'n looks funny b/c we use 
      diploids to refer both to the diploids in a single-deme sim
      or the diploids in a particular deme of a meta-pop sim.

      \note Do not call directly.  Use KTfwd::add_mutation instead.

      \ingroup sugar
    */
    void add_mutation_details( poptype & p,
			       dipvector_t & diploids,
			       const std::size_t mindex,
			       const std::vector<std::size_t> & indlist,
			       const std::vector<short> & clist)
    {
      auto gam_recycling_bin = fwdpp_internal::make_gamete_queue(p.gametes);

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
	      assert(p.gametes[diploids[indlist[i]].first].n); //assert that gamete exists
	      p.gametes[diploids[indlist[i]].first].n--;       //count goes down by 1
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
	      p.gametes[diploids[indlist[i]].first].n++; //update gamete count
	      p.mcounts[mindex]++;
	    }
	  if(clist[i]>0)
	    {
	      assert(p.gametes[diploids[indlist[i]].second].n); //assert that gamete exists
	      p.gametes[diploids[indlist[i]].second].n--;       //count goes down by 1
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
	      p.gametes[diploids[indlist[i]].second].n++; //update gamete count
	      p.mcounts[mindex]++;
	    }
	}
    }

    template<typename poptype,
	     typename dipvector_t>
    /*!
      \brief Overload for adding mutations into single-deme multi-locus populations

      \note Do not call directly.  Use KTfwd::add_mutation instead.

      \ingroup sugar
    */
    void add_mutation_details( poptype & p,
			       dipvector_t & diploids,
			       const std::size_t locus,
			       const std::size_t mindex,
			       const std::vector<std::size_t> & indlist,
			       const std::vector<short> & clist)
    {
      auto gam_recycling_bin = fwdpp_internal::make_gamete_queue(p.gametes);


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
	      assert(p.gametes[diploids[indlist[i]][locus].first].n); //assert that gamete exists
	      p.gametes[diploids[indlist[i]][locus].first].n--;       //count goes down by 1
	      auto n = p.gametes[diploids[indlist[i]][locus].first].mutations;
	      auto s = p.gametes[diploids[indlist[i]][locus].first].smutations;
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
	      diploids[indlist[i]][locus].first = fwdpp_internal::recycle_gamete(p.gametes,
										 gam_recycling_bin,
										 n,s);
	      p.gametes[diploids[indlist[i]][locus].first].n++; //update gamete count
	      p.mcounts[mindex]++;
	    }
	  if(clist[i]>0)
	    {
	      assert(p.gametes[diploids[indlist[i]][locus].second].n); //assert that gamete exists
	      p.gametes[diploids[indlist[i]][locus].second].n--;       //count goes down by 1
	      auto n = p.gametes[diploids[indlist[i]][locus].second].mutations;
	      auto s = p.gametes[diploids[indlist[i]][locus].second].smutations;
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
	      diploids[indlist[i]][locus].second = fwdpp_internal::recycle_gamete(p.gametes,
										  gam_recycling_bin,
										  n,s);
	      p.gametes[diploids[indlist[i]][locus].second].n++; //update gamete count
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

    Note that \a args can take on a few different forms.  First, it can be a raw set of values
    used to construct a new mutation.  Or, it can be an object of correct mutation type.  Or, it can be
    any type from which the correct mutation type can be constructed.  The last two cases require 
    that the mutation type have the appropriate constructors defined.

    See the unit test file unit/test_sugar_add_mutation.cc for example of use.

    \ingroup sugar
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
    auto mindex = sugar::get_mut_index(p.mutations,p.mcounts,new_mutant);
    sugar::add_mutation_details(p,p.diploids,mindex,indlist,clist);
  }

  template<typename metapoptype,
	   class... Args>
  void add_mutation(metapoptype & p,
		    const std::vector<std::size_t> demes,
		    const std::vector<std::vector<std::size_t>> & indlist,
		    const std::vector<std::vector<short>> & clist,
		    Args&&... args)
  /*!
    \brief Add a mutation into a deme from a population at a given frequency in a specific set of
    demes.

    \param p A population object. Meta- or multi-locus.
    \parm demes Vector of indices of the demes in which to add mutation.
    \param indlist A list of indexes of diploids into which to add the new mutations in each deme.
    \param clist A list of gametes. See below.
    \param args Values required to cosnstruct a new mutation.  See below.

    \return Nothing (void)

    Some notes:

    clist.size() must equal indlist.size() must equal demes.size()

    Further, each element of clist must equal each element of indlist in size.

    Values in \a clist must be 0, 1, or 2. These values mean to add the mutation to the first,
    second, or both gametes, resepectively, of each diploid in \a indlist.

    Note that \a args can take on a few different forms.  First, it can be a raw set of values
    used to construct a new mutation.  Or, it can be an object of correct mutation type.  Or, it can be
    any type from which the correct mutation type can be constructed.  The last two cases require 
    that the mutation type have the appropriate constructors defined.

    See the unit test file unit/test_sugar_add_mutation.cc for example of use.

    \ingroup sugar
  */
  {
    for( const auto & c : clist )
      {
	for( auto ci : c )
	  if(ci<0||ci>2) throw std::runtime_error("clist contains elements < 0 and/or > 2");
      }
    if(indlist.size()!=clist.size()) throw std::runtime_error("indlist and clist must be same length");
    for(std::size_t i=0;i<indlist.size();++i)
      {
	if(indlist[i].size()!=clist[i].size()) throw std::runtime_error("indlist[i] must equal clist[i] for all i");
      }

    for(std::size_t i=0;i<demes.size();++i)
      {
	if(i>p.diploids.size()) throw std::out_of_range("deme index out of range");
	auto d = demes[i];
	for( const auto di : indlist[i] )
	  {
	    if(di>p.diploids[d].size()) throw std::out_of_range("individual index out of range");
	  }
      }
    //create a new mutation
    typename metapoptype::mcont_t::value_type new_mutant(std::forward<Args>(args)...);
    auto mindex = sugar::get_mut_index(p.mutations,p.mcounts,new_mutant);
    for( std::size_t i = 0 ; i < demes.size() ; ++i )
      {
	sugar::add_mutation_details(p,p.diploids[demes[i]],mindex,indlist[i],clist[i]);
      }
  }
  
  template<typename multiloc_poptype,
	   class... Args>
  void add_mutation(multiloc_poptype & p,
		    const std::size_t locus,
		    const std::vector<std::size_t> & indlist,
		    const std::vector<short> & clist,
		    Args&&... args)
  /*!
    \brief Add a mutation into a population at a given frequency at in a specific locus.

    \param p A population object. Meta- or multi-locus.
    \parm locus Index of the locus in which to add mutation.
    \param indlist A list of indexes of diploids into which to add the new mutations.
    \param clist A list of gametes. See below.
    \param args Values required to cosnstruct a new mutation.  See below.

    \return Nothing (void)

    Some notes:

    clist.size() must equal indlist.size()

    Values in \a clist must be 0, 1, or 2. These values mean to add the mutation to the first,
    second, or both gametes, resepectively, of each diploid in \a indlist.

    Note that \a args can take on a few different forms.  First, it can be a raw set of values
    used to construct a new mutation.  Or, it can be an object of correct mutation type.  Or, it can be
    any type from which the correct mutation type can be constructed.  The last two cases require 
    that the mutation type have the appropriate constructors defined.

    See the unit test file unit/test_sugar_add_mutation.cc for example of use.

    \ingroup sugar
  */
  {
    for( const auto & c : clist )
      {
	if(c<0||c>2) throw std::runtime_error("clist contains elements < 0 and/or > 2");
      }
    if(indlist.size()!=clist.size()) throw std::runtime_error("indlist and clist must be same length");

    if(locus >= p.diploids[0].size()) throw std::out_of_range("locus index out of range");
    for(const auto & i : indlist)
      {
	if(i>=p.diploids.size()) throw std::out_of_range("indlist contains elements > p.diploids.size()");
      }
    //create a new mutation
    typename multiloc_poptype::mcont_t::value_type new_mutant(std::forward<Args>(args)...);
    auto mindex = sugar::get_mut_index(p.mutations,p.mcounts,new_mutant);
    sugar::add_mutation_details(p,p.diploids,locus,mindex,indlist,clist);
  }
}

#endif
