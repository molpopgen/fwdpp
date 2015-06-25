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
  template< template<typename,typename> class diploid_geno_t,
	    typename gamete_type,
	    typename gamete_list_type_allocator,
	    template<typename,typename> class gamete_list_type,
	    typename mutation_type,
	    typename mutation_list_type_allocator,
	    template<typename,typename> class mutation_list_type,
	    typename vector_type_allocator,
	    template<typename,typename> class diploid_vector_type,
	    typename mutation_writer_type,
	    typename ostreamtype,
	    typename diploid_writer_t = diploidIOplaceholder>
  void write_binary_pop ( const gamete_list_type< gamete_type, gamete_list_type_allocator > * gametes,
			  const mutation_list_type< mutation_type, mutation_list_type_allocator > * mutations,
			  const diploid_vector_type< diploid_geno_t<typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator,typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator>,vector_type_allocator > * diploids,
			  const mutation_writer_type & mw,
			  ostreamtype & buffer,
			  const diploid_writer_t & dw = diploid_writer_t() )
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
  
  template< template<typename,typename> class diploid_geno_t,
	    typename gamete_type,
	    typename gamete_list_type_allocator,
	    template<typename,typename> class gamete_list_type,
	    typename mutation_type,
	    typename mutation_list_type_allocator,
	    template<typename,typename> class mutation_list_type,
	    typename vector_type_allocator,
	    template<typename,typename> class diploid_vector_type,
	    typename mutation_reader_type,
	    typename istreamtype,
	    typename diploid_reader_t=diploidIOplaceholder>
  void read_binary_pop (  gamete_list_type< gamete_type, gamete_list_type_allocator > * gametes,
			  mutation_list_type< mutation_type, mutation_list_type_allocator > * mutations,
			  diploid_vector_type< diploid_geno_t<typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator,typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator>,vector_type_allocator > * diploids,
			  const mutation_reader_type & mr,
			  istreamtype & in,
			  const diploid_reader_t & dr = diploid_reader_t())
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
  template< template<typename,typename> class diploid_geno_t,
	    typename gamete_type,
	    typename gamete_list_type_allocator,
	    template<typename,typename> class gamete_list_type,
	    typename mutation_type,
	    typename mutation_list_type_allocator,
	    template<typename,typename> class mutation_list_type,
	    typename vector_type_allocator,
	    template<typename,typename> class diploid_vector_type,
	    typename diploid_vv_type_allocator,
	    template<typename,typename> class diploid_vv_type,
	    typename mutation_writer_type,
	    typename ostreamtype,
	    typename diploid_writer_t= diploidIOplaceholder>
  void  write_binary_pop ( //const mlocus_vector_type< gamete_list_type< gamete_type, gamete_list_type_allocator >, mlocus_vector_type_allocator> * mlocus_gametes,
			  const gamete_list_type<gamete_type, gamete_list_type_allocator > * mlocus_gametes,
			  const mutation_list_type< mutation_type, mutation_list_type_allocator > * mutations,
			  const diploid_vv_type <diploid_vector_type<diploid_geno_t<typename gamete_list_type<gamete_type, gamete_list_type_allocator >::iterator,typename gamete_list_type<gamete_type,gamete_list_type_allocator>::iterator>, vector_type_allocator >,  diploid_vv_type_allocator > * diploids,
			  const mutation_writer_type & mw,
			  ostreamtype & buffer,
			  const diploid_writer_t & dw = diploid_writer_t() )
  {
    unsigned nloci = mlocus_gametes->size();
    buffer.write(reinterpret_cast<char*>(&nloci),sizeof(unsigned));
    //write mutations
    auto mutdata = fwdpp_internal::write_mutations()(mutations,mw,buffer);
    auto gamdata = fwdpp_internal::write_haplotypes()( mlocus_gametes, mutdata.first, mutdata.second, buffer );
    unsigned ndips=diploids->size();
    buffer.write( reinterpret_cast<char*>(&ndips),sizeof(unsigned) );
    using mloc_diploid_geno_t = typename diploid_vv_type < diploid_vector_type< diploid_geno_t<typename  gamete_list_type<gamete_type, gamete_list_type_allocator >::iterator,typename  gamete_list_type<gamete_type, gamete_list_type_allocator >::iterator>, vector_type_allocator >,  diploid_vv_type_allocator >::value_type;
    std::for_each( diploids->cbegin(), diploids->cend(),
    		   [&gamdata,&buffer,&dw]( const mloc_diploid_geno_t & diploid ) {
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
    // //write haplotypes
    // using gmap_t = std::vector< typename gamete_list_type< gamete_type, gamete_list_type_allocator >::const_iterator >;
    // std::vector< std::pair<gmap_t, std::vector<unsigned> > > gamdata_vec;
    // std::for_each( mlocus_gametes->cbegin(), mlocus_gametes->cend(),
    // 		   [&gamdata_vec,&mutdata,&buffer](const gamete_list_type< gamete_type, gamete_list_type_allocator > & gametes ) {
    // 		     gamdata_vec.emplace_back( std::move( fwdpp_internal::write_haplotypes()(&gametes, mutdata.first,mutdata.second,buffer) ) );
    // 		   }
    // 		   );
    // //Write the diploids
    // unsigned ndips = diploids->size();
    // buffer.write( reinterpret_cast<char*>(&ndips),sizeof(unsigned) );
    // using mloc_diploid_geno_t = typename diploid_vv_type < diploid_vector_type< diploid_geno_t, vector_type_allocator >,  diploid_vv_type_allocator >::value_type;
    // std::for_each( diploids->cbegin(), diploids->cend(),
    // 		   [&gamdata_vec,&buffer,&dw]( const mloc_diploid_geno_t & diploid ) {
    // 		     unsigned i = 0;
    // 		     for( auto genotype = diploid.cbegin() ; genotype != diploid.cend(); ++genotype,++i )
    // 		       {
    // 			 unsigned c = gamdata_vec[i].second[ std::vector<unsigned>::size_type(std::find( gamdata_vec[i].first.begin(),gamdata_vec[i].first.end(),genotype->first ) - gamdata_vec[i].first.begin()) ];
    // 			 buffer.write( reinterpret_cast<char*>(&c),sizeof(unsigned) );
    // 			 c = gamdata_vec[i].second[ std::vector<unsigned>::size_type(std::find( gamdata_vec[i].first.begin(),gamdata_vec[i].first.end(),genotype->second ) - gamdata_vec[i].first.begin()) ];
    // 			 buffer.write( reinterpret_cast<char*>(&c),sizeof(unsigned) );
    // 			 dw(genotype,buffer);
    // 		       }
    // 		   } 
    // 		   );
  }
  
  //Multilocus, single-population, istream
  template< template<typename,typename> class diploid_geno_t,
	    typename gamete_type,
	    typename gamete_list_type_allocator,
	    template<typename,typename> class gamete_list_type,
	    //typename mlocus_vector_type_allocator,
	    //template<typename,typename> class mlocus_vector_type,
	    typename mutation_type,
	    typename mutation_list_type_allocator,
	    template<typename,typename> class mutation_list_type,
	    typename vector_type_allocator,
	    template<typename,typename> class diploid_vector_type,
	    typename diploid_vv_type_allocator,
	    template<typename,typename> class diploid_vv_type,
	    typename mutation_reader_type,
	    typename istreamtype,
	    typename diploid_reader_t= diploidIOplaceholder>
  void read_binary_pop ( //mlocus_vector_type< gamete_list_type< gamete_type, gamete_list_type_allocator >, mlocus_vector_type_allocator> * mlocus_gametes,
			gamete_list_type< gamete_type, gamete_list_type_allocator> * mlocus_gametes,
			 mutation_list_type< mutation_type, mutation_list_type_allocator > * mutations,
			diploid_vv_type < diploid_vector_type< diploid_geno_t<typename  gamete_list_type<gamete_type, gamete_list_type_allocator >::iterator,typename  gamete_list_type<gamete_type, gamete_list_type_allocator >::iterator> , vector_type_allocator >, diploid_vv_type_allocator > * diploids,
			 const mutation_reader_type & mr,
			 istreamtype & in,
			const diploid_reader_t & dr = diploid_reader_t() )
  {
    mlocus_gametes->clear();
    mutations->clear();
    diploids->clear();

    unsigned nloci;
    fwdpp_internal::scalar_reader<unsigned>()(in,&nloci);
    //Read the mutations from the buffer
    auto mutdata = fwdpp_internal::read_mutations()( mutations,mr,in);
    auto gam_info_vec = fwdpp_internal::read_haplotypes()(mlocus_gametes,mutdata,in);
    using mloc_diploid_geno_t = typename diploid_vv_type < diploid_vector_type< diploid_geno_t<typename  gamete_list_type<gamete_type, gamete_list_type_allocator >::iterator,typename  gamete_list_type<gamete_type, gamete_list_type_allocator >::iterator> , vector_type_allocator >, diploid_vv_type_allocator >::value_type;
    unsigned ndips;
    fwdpp_internal::scalar_reader<unsigned>()(in,&ndips);
    diploids->resize(ndips, mloc_diploid_geno_t(nloci,diploid_geno_t<typename  gamete_list_type<gamete_type, gamete_list_type_allocator >::iterator,typename  gamete_list_type<gamete_type, gamete_list_type_allocator >::iterator>()) );
    std::for_each( diploids->begin(), diploids->end(),
		   [&gam_info_vec,&in,&dr]( mloc_diploid_geno_t & diploid ) {
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

    // using gam_info_t = std::map<unsigned,typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator>;
    // std::vector< gam_info_t > gam_info_vec;
    // mlocus_gametes->resize(nloci);
    // //Read the haplotypes
    // for( auto l = mlocus_gametes->begin();l != mlocus_gametes->end(); ++l )
    //   {
    // 	gam_info_vec.emplace_back( std::move( fwdpp_internal::read_haplotypes()(&*l,mutdata,in) ) );
    //   }
    // //read diploids
    // using mloc_diploid_geno_t = typename diploid_vv_type < diploid_vector_type< diploid_geno_t , vector_type_allocator >, diploid_vv_type_allocator >::value_type;
    
    // unsigned ndips;
    // fwdpp_internal::scalar_reader<unsigned>()(in,&ndips);
    // diploids->resize(ndips, mloc_diploid_geno_t(nloci,diploid_geno_t()) );
    // std::for_each( diploids->begin(), diploids->end(),
    // 		   [&gam_info_vec,&in,&dr]( mloc_diploid_geno_t & diploid ) {
    // 		     unsigned i = 0;
    // 		     for( auto l = diploid.begin(); l != diploid.end() ; ++l,++i )
    // 		       {
    // 			 unsigned c;
    // 			 fwdpp_internal::scalar_reader<unsigned>()(in,&c);
    // 			 l->first = gam_info_vec[i][c];
    // 			 fwdpp_internal::scalar_reader<unsigned>()(in,&c);
    // 			 l->second = gam_info_vec[i][c];
    // 			 dr(l,in);
    // 		       }
    // 		   }
    // 		   );
  }

  template< typename diploid_geno_t,
	    typename gamete_type,
	    typename gamete_list_type_allocator,
	    template<typename,typename> class gamete_list_type,
	    typename mutation_type,
	    typename mutation_list_type_allocator,
	    template<typename,typename> class mutation_list_type,
	    typename vector_type_allocator,
	    template<typename,typename> class diploid_vector_type,
	    typename diploid_vv_type_allocator,
	    template<typename,typename> class diploid_vv_type,
	    typename mutation_writer_type,
	    typename ostreamtype,
	    typename diploid_writer_t= diploidIOplaceholder>
  void write_binary_metapop (const gamete_list_type< gamete_type, gamete_list_type_allocator > * metapop,
			     const mutation_list_type< mutation_type, mutation_list_type_allocator > * mutations,
			     const diploid_vv_type < diploid_vector_type< diploid_geno_t , vector_type_allocator >, diploid_vv_type_allocator > * diploids,
			     const mutation_writer_type & mw,
			     ostreamtype & buffer,
			     const diploid_writer_t & dw = diploid_writer_t())
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
  
  template< typename diploid_geno_t,
	    typename gamete_type,
  	    typename gamete_list_type_allocator,
  	    template<typename,typename> class gamete_list_type,
  	    typename mutation_type,
  	    typename mutation_list_type_allocator,
  	    template<typename,typename> class mutation_list_type,
  	    typename vector_type_allocator,
  	    template<typename,typename> class diploid_vector_type,
  	    typename diploid_vv_type_allocator,
  	    template<typename,typename> class diploid_vv_type,
  	    typename mutation_reader_type,
  	    typename istreamtype,
	    typename diploid_reader_t= diploidIOplaceholder>
  void read_binary_metapop (gamete_list_type< gamete_type, gamete_list_type_allocator > * metapop,
			    mutation_list_type< mutation_type, mutation_list_type_allocator > * mutations,
			    diploid_vv_type < diploid_vector_type<diploid_geno_t ,vector_type_allocator >, diploid_vv_type_allocator > * diploids,
			    const mutation_reader_type & mr,
			    istreamtype & in,
			    const diploid_reader_t & dr = diploid_reader_t() )
  {
    metapop->clear();
    mutations->clear();
    diploids->clear();
    
    using mlist = mutation_list_type< mutation_type, mutation_list_type_allocator >;
    using mut_info = std::map<unsigned,typename mlist::iterator>;
    
    unsigned NPOP;
    fwdpp_internal::scalar_reader<unsigned>()(in,&NPOP);
    mut_info m = fwdpp_internal::read_mutations()(mutations,mr,in); 
    
    diploids->resize(NPOP);
    auto dip_ptr = diploids->begin();
    using glist = gamete_list_type< gamete_type, gamete_list_type_allocator >;
    using gam_info = std::map<unsigned,typename glist::iterator>;
    gam_info g = fwdpp_internal::read_haplotypes()(metapop,m,in);	
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
