//  -*- C++ -*- 
#ifndef __FWDPP_IO_TCC__
#define __FWDPP_IO_TCC__

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
  template< typename gcont_t,
	    typename mcont_t,
	    typename dipvector_t,
	    typename mutation_writer_type,
	    typename ostreamtype,
	    typename diploid_writer_t>
  void write_binary_pop ( const gcont_t & gametes,
			  const mcont_t & mutations,
			  const dipvector_t & diploids,
			  const mutation_writer_type & mw,
			  ostreamtype & buffer,
			  const diploid_writer_t & dw )
  {
    fwdpp_internal::write_mutations()( mutations,mw,buffer ); 
    fwdpp_internal::write_haplotypes()( gametes, buffer );
    std::size_t NDIPS = diploids.size();
    buffer.write( reinterpret_cast<char *>(&NDIPS), sizeof(std::size_t) );
    for(const auto & dip : diploids )
      {
	buffer.write(reinterpret_cast<const char *>(&dip.first),sizeof(std::size_t));
	buffer.write(reinterpret_cast<const char *>(&dip.second),sizeof(std::size_t));
	dw(dip,buffer);
      }
  }
  
  template< typename gcont_t,
	    typename mcont_t,
	    typename dipvector_t,
	    typename mutation_reader_type,
	    typename istreamtype,
	    typename diploid_reader_t>
  void read_binary_pop (  gcont_t & gametes,
			  mcont_t & mutations,
			  dipvector_t & diploids,
			  const mutation_reader_type & mr,
			  istreamtype & in,
			  const diploid_reader_t & dr)
  {
    gametes.clear();
    mutations.clear();
    diploids.clear();
    fwdpp_internal::read_mutations()(mutations,mr,in); 
    fwdpp_internal::read_haplotypes()(gametes,in);
    std::size_t NDIPS,c;
    fwdpp_internal::scalar_reader<decltype(NDIPS)>()(in,&NDIPS);
    diploids.resize(NDIPS);
    for( auto & dip : diploids )
      {
	fwdpp_internal::scalar_reader<decltype(c)>()(in,&c);
	dip.first = c;
	fwdpp_internal::scalar_reader<decltype(c)>()(in,&c);
	dip.second = c;
	dr(dip,in);
      }
  }

  //multi-locus, single pop, ostream
  template< typename gcont_t,
	    typename mcont_t,
	    typename dipvector_t,
	    typename mutation_writer_type,
	    typename ostreamtype,
	    typename diploid_writer_t>
  void  write_binary_pop_mloc ( const gcont_t & mlocus_gametes,
				const mcont_t & mutations,
				const dipvector_t & diploids,
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
  template< typename gcont_t,
	    typename mcont_t,
	    typename dipvector_t,
	    typename mutation_reader_type,
	    typename istreamtype,
	    typename diploid_reader_t>
  void read_binary_pop_mloc (gcont_t & mlocus_gametes,
			     mcont_t & mutations,
			     dipvector_t & diploids,
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

  template< typename gcont_t,
	    typename mcont_t,
	    typename dipvector_t,
	    typename mutation_writer_type,
	    typename ostreamtype,
	    typename diploid_writer_t>
  void write_binary_metapop (const gcont_t & gametes,
			     const mcont_t & mutations,
			     const dipvector_t & diploids,
			     const mutation_writer_type & mw,
			     ostreamtype & buffer,
			     const diploid_writer_t & dw)
  {
    std::size_t i = unsigned(diploids.size());
    buffer.write( reinterpret_cast<char *>(&i), sizeof(decltype(i)) );
    fwdpp_internal::write_mutations()( mutations,mw,buffer );
    fwdpp_internal::write_haplotypes()( gametes, buffer );
    for(const auto & deme : diploids)
      {
	i = deme.size();
	buffer.write( reinterpret_cast<char *>(&i), sizeof(decltype(i)) );
	for(const auto & dip : deme)
	  {
	    buffer.write(reinterpret_cast<const char*>(&dip.first),sizeof(decltype(dip.first)));
	    buffer.write(reinterpret_cast<const char*>(&dip.second),sizeof(decltype(dip.second)));
	    dw(dip,buffer);
	  }
      }
  }
  
  template< typename gcont_t,
	    typename mcont_t,
	    typename dipvector_t,
  	    typename mutation_reader_type,
  	    typename istreamtype,
	    typename diploid_reader_t>
  void read_binary_metapop (gcont_t & gametes,
			    mcont_t & mutations,
			    dipvector_t & diploids,
			    const mutation_reader_type & mr,
			    istreamtype & in,
			    const diploid_reader_t & dr)
  {
    gametes.clear();
    mutations.clear();
    diploids.clear();
    
    std::size_t i;
    fwdpp_internal::scalar_reader<decltype(i)>()(in,&i);
    diploids.resize(i);
    fwdpp_internal::read_mutations()(mutations,mr,in); 
    fwdpp_internal::read_haplotypes()(gametes,in);	
    for( auto & deme : diploids)
      {
	fwdpp_internal::scalar_reader<decltype(i)>()(in,&i);
	if(i)
	  {
	    deme.resize(i);
	    for(auto & dip : deme)
	      {
		fwdpp_internal::scalar_reader<decltype(dip.first)>()(in,&dip.first);
		fwdpp_internal::scalar_reader<decltype(dip.second)>()(in,&dip.second);
		dr(dip,in);
	      }
	  }
      }
  }

}//ns KTfwd

#endif
