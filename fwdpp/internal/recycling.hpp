#ifndef FWDPP_INTERNAL_RECYCLING
#define FWDPP_INTERNAL_RECYCLING

#include <queue>

namespace KTfwd
{
  namespace fwdpp_internal
  {
    template<class T> using recycling_bin_t = std::queue<T>;

    template<typename mcount_vec>
    recycling_bin_t<std::size_t> make_mut_queue( const mcount_vec & mcounts )
    {
      recycling_bin_t<std::size_t> rv;
      for( std::size_t i=0;i<mcounts.size();++i)
	{
	  if(!mcounts[i])rv.push(i);
	}
      return rv;
    }

    template<typename gvec_t>
    recycling_bin_t<std::size_t> make_gamete_queue( const gvec_t & gametes )
    {
      recycling_bin_t<std::size_t> rv;
      for(std::size_t i=0;i<gametes.size();++i)
	{
	  if(!gametes[i].n)rv.push(i);
	}
      return rv;
    }

    template<typename gcont_t,
	     typename mcont_t,
	     typename queue_t,
	     typename lookup_type>
    inline std::size_t recycle_gamete(	gcont_t & gametes,
					const mcont_t & mutations,
					queue_t & gamete_recycling_bin, lookup_type & gamete_lookup,
					std::vector<std::size_t> & neutral,
					std::vector<std::size_t> & selected )
    {
      //Try to recycle
      if( ! gamete_recycling_bin.empty() )
	{
	  auto idx = gamete_recycling_bin.front();
	  gamete_recycling_bin.pop();
	  assert(!gametes[idx].n);
	  gametes[idx].n=0u;
	  gametes[idx].mutations.swap(neutral);
	  gametes[idx].smutations.swap(selected);
	  gamete_lookup.update(idx,gametes,mutations);
	  return idx;
	}
      //gametes.emplace_back(0u,std::move(neutral),std::move(selected));
      gametes.emplace_back(0u,neutral,selected);
      gamete_lookup.update(gametes.size()-1,gametes,mutations);
      assert( gametes[gametes.size()-1].mutations==neutral );
      assert( gametes[gametes.size()-1].smutations==selected );
      return (gametes.size()-1);
    }

    /*!
      \brief Helper function for mutation policies

      This function minimizes code duplication when writing mutation models.  It abstracts
      the operations needed to recycle an extinct mutation.

      \param mutation_recycling_bin  A FIFO queue of iterators pointing to extinct mutations.
      \param mutations A list of mutation objects
      \param args Parameter pack to be passed to constructor of an mlist_t::value_type
     */
    template<typename queue_t,
	     typename mlist_t,
	     class... Args >
    typename std::size_t recycle_mutation_helper( queue_t & mutation_recycling_bin,
						  mlist_t & mutations,
						  Args&&... args )
    {
      if(!mutation_recycling_bin.empty())
	{
	  auto rv = mutation_recycling_bin.front();
	  mutation_recycling_bin.pop();
	  mutations[rv]=typename mlist_t::value_type(args...);
	  return rv;
	}
      mutations.emplace_back(std::forward<Args>(args)...);
      return mutations.size()-1;
    }
  }
}

#endif
