#ifndef __FWDPP_INTERNAL_MS_SAMPLING_HPP__
#define __FWDPP_INTERNAL_MS_SAMPLING_HPP__

namespace KTfwd {
  namespace fwdpp_internal {
    template< typename mcont_t,
	      typename pos_finder>
    void update_sample_block(std::vector< std::pair<double,std::string> > & block,
			     const mcont_t & mutations,
			     const unsigned & i,
			     const unsigned & n,
			     const pos_finder & pf,
			     const unsigned & offset = 0)
    {
      for( auto mptr = mutations.cbegin() ; mptr != mutations.cend() ; ++mptr )
	{
	  double mutpos = (*mptr)->pos;
	  auto itr = std::find_if(block.begin(),block.end(),std::bind(pf,std::placeholders::_1,mutpos));
	  if( itr == block.end() )
	    {
	      block.push_back( std::make_pair(mutpos,std::string(n,'0')) );
	      block[block.size()-1].second[2*i+offset] = '1';
	    }
	  else
	    {
	      itr->second[2*i+offset]='1';
	    }
	}
    }
  }
}

#endif
