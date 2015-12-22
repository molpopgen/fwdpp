//  -*- C++ -*- 
#ifndef __FWDPP_IO_TCC__
#define __FWDPP_IO_TCC__

#include <map>
#include <vector>
#include <algorithm>
#include <cassert>
#include <type_traits>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/internal/IOhelp.hpp>

namespace KTfwd
{  
  //Binary I/O for individual-based simulation

  //Single-locus sims, single pop
  template< typename glist_t,
	    typename mlist_t,
	    typename dipvector_t,
	    typename mutation_writer_type,
	    typename ostreamtype,
	    typename diploid_writer_t>
  void write_binary_pop ( const glist_t * gametes,
			  const mlist_t * mutations,
			  const dipvector_t * diploids,
			  const mutation_writer_type & mw,
			  ostreamtype & buffer,
			  const diploid_writer_t & dw )
  {
    auto mutdata = fwdpp_internal::write_mutations()( mutations,mw,buffer ); 
    auto gamdata = fwdpp_internal::write_haplotypes()( gametes, mutdata.first, mutdata.second, buffer );
    unsigned NDIPS = unsigned(diploids->size());
    buffer.write( reinterpret_cast<char *>(&NDIPS), sizeof(unsigned) );
    for( auto dip = diploids->cbegin() ; dip != diploids->cend() ; ++dip )
      {
	unsigned c = gamdata.second[ std::vector<unsigned>::size_type(std::find( gamdata.first.begin(),gamdata.first.end(),dip->first ) - gamdata.first.begin()) ];
	buffer.write( reinterpret_cast<char *>(&c), sizeof(unsigned) );
	c = gamdata.second[ std::vector<unsigned>::size_type(std::find( gamdata.first.begin(),gamdata.first.end(),dip->second ) - gamdata.first.begin()) ];
	buffer.write( reinterpret_cast<char *>(&c), sizeof(unsigned) );
	dw(dip,buffer);
      }
  }
  
  template< typename glist_t,
	    typename mlist_t,
	    typename dipvector_t,
	    typename mutation_reader_type,
	    typename istreamtype,
	    typename diploid_reader_t>
  void read_binary_pop (  glist_t * gametes,
			  mlist_t * mutations,
			  dipvector_t * diploids,
			  const mutation_reader_type & mr,
			  istreamtype & in,
			  const diploid_reader_t & dr)
  {
    gametes->clear();
    mutations->clear();
    diploids->clear();
    auto m = fwdpp_internal::read_mutations()(mutations,mr,in); 
    auto g = fwdpp_internal::read_haplotypes()(gametes,m,in);
    unsigned NDIPS,c;
    fwdpp_internal::scalar_reader<unsigned>()(in,&NDIPS);
    diploids->resize(NDIPS);
    for( auto dp = diploids->begin() ; dp != diploids->end() ; ++dp )
      {
	fwdpp_internal::scalar_reader<unsigned>()(in,&c);
	dp->first = g[c];
	fwdpp_internal::scalar_reader<unsigned>()(in,&c);
	dp->second = g[c];
	dr(dp,in);
      }
  }

  //multi-locus, single pop, ostream
  template< typename glist_t,
	    typename mlist_t,
	    typename dipvector_t,
	    typename mutation_writer_type,
	    typename ostreamtype,
	    typename diploid_writer_t>
  void  write_binary_pop_mloc ( const glist_t * mlocus_gametes,
				const mlist_t * mutations,
				const dipvector_t * diploids,
				const mutation_writer_type & mw,
				ostreamtype & buffer,
				const diploid_writer_t & dw )
  {
    unsigned nloci = unsigned(diploids->begin()->size());
    buffer.write(reinterpret_cast<char*>(&nloci),sizeof(unsigned));
    //write mutations
    auto mutdata = fwdpp_internal::write_mutations()(mutations,mw,buffer);
    auto gamdata = fwdpp_internal::write_haplotypes()( mlocus_gametes, mutdata.first, mutdata.second, buffer );
    unsigned ndips=unsigned(diploids->size());
    buffer.write( reinterpret_cast<char*>(&ndips),sizeof(unsigned) );
    std::for_each( diploids->cbegin(), diploids->cend(),
    		   [&gamdata,&buffer,&dw]( const typename dipvector_t::value_type & diploid ) {
    		     unsigned i = 0;
    		     for( auto genotype = diploid.cbegin() ; genotype != diploid.cend(); ++genotype,++i )
    		       {
    			 unsigned c = gamdata.second[ std::vector<unsigned>::size_type(std::find( gamdata.first.begin(),gamdata.first.end(),genotype->first ) - gamdata.first.begin()) ];
    			 buffer.write( reinterpret_cast<char*>(&c),sizeof(unsigned) );
    			 c = gamdata.second[ std::vector<unsigned>::size_type(std::find( gamdata.first.begin(),gamdata.first.end(),genotype->second ) - gamdata.first.begin()) ];
    			 buffer.write( reinterpret_cast<char*>(&c),sizeof(unsigned) );
    			 dw(genotype,buffer);
    		       }
    		   } 
    		   );
  }
  
  //Multilocus, single-population, istream
  template< typename glist_t,
	    typename mlist_t,
	    typename dipvector_t,
	    typename mutation_reader_type,
	    typename istreamtype,
	    typename diploid_reader_t>
  void read_binary_pop_mloc (glist_t * mlocus_gametes,
			     mlist_t * mutations,
			     dipvector_t * diploids,
			     const mutation_reader_type & mr,
			     istreamtype & in,
			     const diploid_reader_t & dr)
  {
    mlocus_gametes->clear();
    mutations->clear();
    diploids->clear();

    unsigned nloci;
    fwdpp_internal::scalar_reader<unsigned>()(in,&nloci);
    //Read the mutations from the buffer
    auto mutdata = fwdpp_internal::read_mutations()( mutations,mr,in);
    auto gam_info_vec = fwdpp_internal::read_haplotypes()(mlocus_gametes,mutdata,in);
    unsigned ndips;
    fwdpp_internal::scalar_reader<unsigned>()(in,&ndips);
    diploids->resize(ndips, typename dipvector_t::value_type(nloci) ); 
    std::for_each( diploids->begin(), diploids->end(),
		   [&gam_info_vec,&in,&dr]( typename dipvector_t::value_type  & diploid ) {
		     unsigned i = 0;
     		     for( auto l = diploid.begin(); l != diploid.end() ; ++l,++i )
     		       {
     			 unsigned c;
     			 fwdpp_internal::scalar_reader<unsigned>()(in,&c);
     			 l->first = gam_info_vec[c];
     			 fwdpp_internal::scalar_reader<unsigned>()(in,&c);
     			 l->second = gam_info_vec[c];
     			 dr(l,in);
     		       }
     		   }
     		   );
  }

  template< typename glist_t,
	    typename mlist_t,
	    typename dipvector_t,
	    typename mutation_writer_type,
	    typename ostreamtype,
	    typename diploid_writer_t>
  void write_binary_metapop (const glist_t * metapop,
			     const mlist_t * mutations,
			     const dipvector_t * diploids,
			     const mutation_writer_type & mw,
			     ostreamtype & buffer,
			     const diploid_writer_t & dw)
  {
    unsigned NPOP = unsigned(diploids->size());
    buffer.write( reinterpret_cast<char *>(&NPOP), sizeof(unsigned) );
    auto mutdata = fwdpp_internal::write_mutations()( mutations,mw,buffer );
    auto gamdata = fwdpp_internal::write_haplotypes()( metapop, mutdata.first, mutdata.second, buffer );
    auto dip_ptr = diploids->begin();
    unsigned NDIPS,c;
    for( unsigned pop = 0 ; pop < NPOP ; ++pop,++dip_ptr )
      {
	NDIPS = unsigned(dip_ptr->size());
      	buffer.write( reinterpret_cast<char *>(&NDIPS),sizeof(unsigned) );
	for( auto dip = dip_ptr->cbegin() ; dip != dip_ptr->cend() ; ++dip )
	  {
	    c = gamdata.second[ std::vector<unsigned>::size_type(std::find( gamdata.first.begin(), gamdata.first.end(), dip->first ) - gamdata.first.begin()) ];
	    buffer.write( reinterpret_cast<char *>(&c),sizeof(unsigned) );
	    c = gamdata.second[ std::vector<unsigned>::size_type(std::find( gamdata.first.begin(), gamdata.first.end(), dip->second ) - gamdata.first.begin()) ];
	    buffer.write( reinterpret_cast<char *>(&c),sizeof(unsigned) );
	    dw(dip,buffer);
	  }
      }
  }
  
  template< typename glist_t,
	    typename mlist_t,
	    typename dipvector_t,
  	    typename mutation_reader_type,
  	    typename istreamtype,
	    typename diploid_reader_t>
  void read_binary_metapop (glist_t * metapop,
			    mlist_t * mutations,
			    dipvector_t * diploids,
			    const mutation_reader_type & mr,
			    istreamtype & in,
			    const diploid_reader_t & dr)
  {
    metapop->clear();
    mutations->clear();
    diploids->clear();
    
    unsigned NPOP;
    fwdpp_internal::scalar_reader<unsigned>()(in,&NPOP);
    auto m = fwdpp_internal::read_mutations()(mutations,mr,in); 
    
    diploids->resize(NPOP);
    auto dip_ptr = diploids->begin();
    auto g = fwdpp_internal::read_haplotypes()(metapop,m,in);	
    unsigned NDIPS,c;
    for( unsigned pop=0 ; pop < NPOP ; ++pop,++dip_ptr )
      {
	fwdpp_internal::scalar_reader<unsigned>()(in,&NDIPS);
      	dip_ptr->resize(NDIPS);
	for( auto dip = dip_ptr->begin() ; dip != dip_ptr->end() ; ++dip )
	  {
	    fwdpp_internal::scalar_reader<unsigned>()(in,&c);
	    dip->first = g[c];
	    fwdpp_internal::scalar_reader<unsigned>()(in,&c);
	    dip->second = g[c];
	    dr(dip,in);
	  }
      }
  }

}//ns KTfwd

#endif
