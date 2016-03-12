#ifndef __FWDPP_INTERNAL_IOHELP_HPP__
#define __FWDPP_INTERNAL_IOHELP_HPP__

/*
  Mechanics of data serialization
  The various write_binary_pop and read_binary pop
  functions rely on these implementations
*/
#include <zlib.h>

namespace KTfwd {
  namespace fwdpp_internal {

    template<typename T>
    struct scalar_reader
    {
      template<typename streamtype>
      inline void operator()( streamtype & i, T * __t ) const
      {
	i.read( reinterpret_cast<char*>(__t), sizeof(T) );
      }
      inline void operator()( gzFile & gzin, T * __t ) const
      {
	gzread(gzin,__t,sizeof(T));
      }
    };
    
    struct write_mutations
    {
      template< typename mutation_type,
		typename list_type_allocator,
		template<typename,typename> class list_type,
		typename mutation_writer_type,
		typename ostreamtype>
      void operator()( const list_type< mutation_type, list_type_allocator > & mutations,
		       const mutation_writer_type & mw,
		       ostreamtype & buffer) const
      {
	std::size_t MUTNO = mutations.size();
	buffer.write( reinterpret_cast<char *>(&MUTNO), sizeof(std::size_t) );
	//write the mutation data to the buffer
	for(const auto & m : mutations ) mw(m,buffer);
      }
    };

    struct write_haplotypes
    {
      template<
	typename gamete_type,
	typename... gamete_cont_t_details,
	template<typename,typename... > class gamete_cont_t,
	typename ostreamtype>
      void operator()( const gamete_cont_t< gamete_type, gamete_cont_t_details... > & gametes,
		       ostreamtype & buffer) const
      {
	std::size_t N = gametes.size();
	buffer.write( reinterpret_cast< char * >(&N), sizeof(std::size_t) );
	for( const auto & g : gametes )
	  {
	    buffer.write(reinterpret_cast<const char *>(&g.n),sizeof(decltype(g.n)));
	    std::size_t nm = g.mutations.size();
	    buffer.write(reinterpret_cast<const char *>(&nm),sizeof(std::size_t));
	    if(nm)
	      {
		buffer.write(reinterpret_cast<const char *>(&g.mutations[0]),
			     nm*sizeof(typename gamete_type::mutation_container::value_type));
	      }
	    nm = g.smutations.size();
	    buffer.write(reinterpret_cast<const char *>(&nm),sizeof(std::size_t));
	    if(nm)
	      {
		buffer.write(reinterpret_cast<const char *>(&g.smutations[0]),
			     nm*sizeof(typename gamete_type::mutation_container::value_type));
	      }

	  }
      }
    };

    struct read_mutations
    {
      template< typename mutation_type,
		typename list_type_allocator,
		template<typename,typename> class list_type,
		typename mutation_reader,
		typename istreamtype >
      void operator()(list_type< mutation_type, list_type_allocator > & mutations,
		      const mutation_reader & mr,
		      istreamtype & in) const
      {
	std::size_t NMUTS;
	in.read( reinterpret_cast<char *>(&NMUTS), sizeof(decltype(NMUTS)) );
	for(uint_t i = 0 ; i < NMUTS ; ++i)
	  {
	    mutations.emplace_back(mr(in));
	  }
      }
    
      template< typename mutation_type,
		typename list_type_allocator,
		template<typename,typename> class list_type,
		typename mutation_reader>
      void
      operator()(list_type< mutation_type, list_type_allocator > & mutations,
		 const mutation_reader & mr,
		 gzFile gzin) const
      {
	std::size_t NMUTS;
	gzread( gzin, &NMUTS, sizeof(decltype(NMUTS)) );
	for(uint_t i = 0 ; i < NMUTS ; ++i)
	  {
	    mutations.emplace_back( mr(gzin) );
	  }
      }
    };
  
    struct read_haplotypes
    {
      template< typename gamete_type,
		typename list_type_allocator,
		template<typename,typename> class list_type,
		typename istreamtype >
      void  operator()(list_type< gamete_type, list_type_allocator > & gametes,
		       istreamtype & in) const
      {
	std::size_t NHAPS;
	in.read( reinterpret_cast<char *>(&NHAPS),sizeof(decltype(NHAPS)) );
	uint_t N;
	std::size_t nm;
	for( uint_t i = 0 ; i < NHAPS ; ++i )
	  {
	    in.read(reinterpret_cast<char *>(&N),sizeof(decltype(N)));
	    gamete_type g(N);
	    in.read(reinterpret_cast<char *>(&nm),sizeof(decltype(nm)));
	    if(nm)
	      {
		g.mutations.resize(nm);
		in.read(reinterpret_cast<char *>(&g.mutations[0]),nm*sizeof(typename gamete_type::mutation_container::value_type));
	      }
	    in.read(reinterpret_cast<char *>(&nm),sizeof(decltype(nm)));
	    if(nm)
	      {
		g.smutations.resize(nm);
		in.read(reinterpret_cast<char *>(&g.smutations[0]),nm*sizeof(typename gamete_type::mutation_container::value_type));
	      }
	    gametes.emplace_back(std::move(g));
	  }
      }

      template< typename gamete_type,
		typename list_type_allocator,
		template<typename,typename> class list_type>
      void operator()(list_type< gamete_type, list_type_allocator > & gametes,
		      gzFile gzin) const
      {
	std::size_t NHAPS;
	gzread( gzin,&NHAPS,sizeof(decltype(NHAPS)) );
	uint_t N;
	std::size_t nm;
	for( uint_t i = 0 ; i < NHAPS ; ++i )
	  {
	    gzread( gzin,&N,sizeof(decltype(N)) );
	    gamete_type g(N);
	    gzread( gzin,&nm,sizeof(decltype(nm)) );
	    if(nm)
	      {
		g.mutations.resize(nm);
		gzread( gzin,&g.mutations[0],nm*sizeof(typename gamete_type::mutation_container::value_type));
	      }
	    gzread( gzin,&nm,sizeof(decltype(nm)) );
	    if(nm)
	      {
		g.smutations.resize(nm);
		gzread( gzin,&g.smutations[0],nm*sizeof(typename gamete_type::mutation_container::value_type));
	      }
	    gametes.emplace_back(std::move(g));
	  }
      }
    };
  }
}

#endif
