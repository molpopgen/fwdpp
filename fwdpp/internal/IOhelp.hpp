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

    struct standard_diploid_writer
    {
      using result_type = void;
      template< typename iterator >
      inline result_type operator()( iterator i ) const
      {
	//Does nothing!
      }
    };
    
    struct write_mutations
    {
      template< typename mutation_type,
		typename list_type_allocator,
		template<typename,typename> class list_type,
		typename mutation_writer_type,
		typename ostreamtype>
      std::pair< std::vector< typename list_type< mutation_type, list_type_allocator >::const_iterator >,
		 std::vector<unsigned> >
      operator()( const list_type< mutation_type, list_type_allocator > * mutations,
		  const mutation_writer_type & mw,
		  ostreamtype & buffer) const
      {
	using mlist_iterator = typename list_type< mutation_type, list_type_allocator >::const_iterator;
	//initiate a list with the mutation information by iterator
	using maptype = std::vector< mlist_iterator >;
	std::vector<unsigned> indexes;
	maptype mut_info;
		
	unsigned MUTNO=0;
	for( mlist_iterator mtr = mutations->begin() ; mtr != mutations->end() ; ++mtr )
	  {
	    mut_info.push_back( mtr );
	    indexes.push_back(MUTNO++);
	  }
		
	buffer.write( reinterpret_cast<char *>(&MUTNO), sizeof(unsigned) );
		
	//write the mutation data to the buffer
	for( typename maptype::const_iterator itr = mut_info.begin() ;
	     itr != mut_info.end() ; ++itr )
	  {
	    unsigned ID = indexes[std::vector<unsigned>::size_type(itr-mut_info.begin())];
	    buffer.write( reinterpret_cast< char * >(&ID), sizeof(unsigned) );
	    mw( **itr,buffer );
	  }
	return std::make_pair( mut_info, indexes );
      }
    };

    struct write_haplotypes
    {
      template<
	typename gamete_type,
	typename... gamete_cont_t_details,
	template<typename,typename... > class gamete_cont_t,
	typename ostreamtype>
      std::pair< std::vector< typename gamete_cont_t< gamete_type, gamete_cont_t_details... >::const_iterator >,
		 std::vector<unsigned> >
      operator()( const gamete_cont_t< gamete_type, gamete_cont_t_details... > * gametes,
		  const std::vector< typename gamete_type::mutation_list_type::const_iterator > & mut_info,
		  const std::vector< unsigned > & indexes,
		  ostreamtype & buffer) const
      {
	using glist_iterator = typename gamete_cont_t< gamete_type, gamete_cont_t_details...>::const_iterator;
	unsigned N = unsigned(gametes->size());
	buffer.write( reinterpret_cast< char * >(&N), sizeof(unsigned) );
	std::vector< glist_iterator > gam_info;
	std::vector<unsigned> gam_indexes;
	unsigned index = 0;
	for( glist_iterator gptr = gametes->begin() ; gptr != gametes->end() ; ++gptr,++index )
	  {
	    gam_info.push_back( gptr );
	    gam_indexes.push_back(index);
	    buffer.write( reinterpret_cast< char * >(&index),sizeof(unsigned) );
	    N = gptr->n;
	    buffer.write( reinterpret_cast< char * >(&N),sizeof(unsigned) );
	    N = unsigned(gptr->mutations.size());
	    buffer.write( reinterpret_cast<char *>(&N), sizeof(unsigned) );
	    for( unsigned i = 0 ; i < N ; ++i )
	      {
		assert( std::find(mut_info.begin(),mut_info.end(),(gptr->mutations[i])) != mut_info.end() );
		unsigned INDEX = indexes[ std::vector<unsigned>::size_type(std::find(mut_info.begin(),mut_info.end(),(gptr->mutations[i])) - mut_info.begin()) ];
		buffer.write( reinterpret_cast< char * >(&INDEX), sizeof(unsigned) );
	      }
	    N = unsigned(gptr->smutations.size());
	    buffer.write( reinterpret_cast<char *>(&N), sizeof(unsigned) );
	    for( unsigned i = 0 ; i < N ; ++i )
	      {
		assert( std::find(mut_info.begin(),mut_info.end(),(gptr->smutations[i])) != mut_info.end() );
		unsigned INDEX = indexes[ std::vector<unsigned>::size_type(std::find(mut_info.begin(),mut_info.end(),(gptr->smutations[i])) - mut_info.begin()) ];
		buffer.write( reinterpret_cast< char * >(&INDEX), sizeof(unsigned) );
	      }	  
	  }
	return std::make_pair( gam_info, gam_indexes );
      }
    };

    struct read_mutations
    {
      template< typename mutation_type,
		typename list_type_allocator,
		template<typename,typename> class list_type,
		typename mutation_reader,
		typename istreamtype >
      std::map< unsigned, typename list_type<mutation_type,list_type_allocator>::iterator >
      operator()(list_type< mutation_type, list_type_allocator > * mutations,
		 const mutation_reader & mr,
		 istreamtype & in) const
      {
	using mlist = list_type< mutation_type, list_type_allocator >;
	using mut_info = std::map<unsigned,typename mlist::iterator>;
      
	mut_info m;
	unsigned NMUTS;
	in.read( reinterpret_cast<char *>(&NMUTS), sizeof(unsigned) );
	for(unsigned i = 0 ; i < NMUTS ; ++i)
	  {
	    unsigned ID;
	    in.read( reinterpret_cast<char *>(&ID), sizeof(unsigned) );
	    m[ID] = mutations->insert(mutations->end(),mr(in));
	  }
	return m;
      }
    
      template< typename mutation_type,
		typename list_type_allocator,
		template<typename,typename> class list_type,
		typename mutation_reader>
      std::map< unsigned, typename list_type<mutation_type,list_type_allocator>::iterator >
      operator()(list_type< mutation_type, list_type_allocator > * mutations,
		 const mutation_reader & mr,
		 gzFile gzin) const
      {
	using mlist = list_type< mutation_type, list_type_allocator >;
	using mut_info = std::map<unsigned,typename mlist::iterator>;
      
	mut_info m;
	unsigned NMUTS;
	gzread( gzin, &NMUTS, sizeof(unsigned) );
	for(unsigned i = 0 ; i < NMUTS ; ++i)
	  {
	    unsigned ID;
	    gzread( gzin, &ID, sizeof(unsigned) );
	    m[ID] = mutations->insert(mutations->end(),mr(gzin));
	  }
	return m;
      }
    };
  
    struct read_haplotypes
    {
      template< typename gamete_type,
		typename list_type_allocator,
		template<typename,typename> class list_type,
		typename istreamtype >
      std::map< unsigned,
		typename list_type< gamete_type, list_type_allocator >::iterator >
      operator()(list_type< gamete_type, list_type_allocator > * gametes,
		 std::map<unsigned, typename gamete_type::mutation_list_type_iterator> & m,
		 istreamtype & in) const
      {
	std::map< unsigned,
		  typename list_type< gamete_type, list_type_allocator >::iterator > rv;
	unsigned NHAPS;
      
	in.read( reinterpret_cast<char *>(&NHAPS),sizeof(unsigned) );
      
	unsigned INDEX,N,ID_J;
	for( unsigned i = 0 ; i < NHAPS ; ++i )
	  {
	    in.read( reinterpret_cast< char * >(&INDEX), sizeof(unsigned) );
	    assert( rv.find(INDEX) == rv.end() );
	    in.read( reinterpret_cast< char * >(&N), sizeof(unsigned) );
	    gamete_type g(N);
	    in.read( reinterpret_cast< char * >(&N), sizeof(unsigned) );
	    for( unsigned j = 0 ; j < N ; ++j )
	      {
		in.read( reinterpret_cast< char * >(&ID_J), sizeof(unsigned) );
		g.mutations.push_back(m[ID_J]);
	      }
	    in.read( reinterpret_cast< char * >(&N), sizeof(unsigned) );
	    for( unsigned j = 0 ; j < N ; ++j )
	      {
		in.read( reinterpret_cast< char * >(&ID_J), sizeof(unsigned) );
		g.smutations.push_back(m[ID_J]);
	      }
	    rv[INDEX] = gametes->insert(gametes->end(),g);
	  }
	return rv;
      }

      template< typename gamete_type,
		typename list_type_allocator,
		template<typename,typename> class list_type>
      std::map< unsigned,
		typename list_type< gamete_type, list_type_allocator >::iterator >
      operator()(list_type< gamete_type, list_type_allocator > * gametes,
		 std::map<unsigned, typename gamete_type::mutation_list_type_iterator> & m,
		 gzFile gzin) const
      {
	std::map< unsigned,
		  typename list_type< gamete_type, list_type_allocator >::iterator > rv;
	unsigned NHAPS;
      
	gzread( gzin,&NHAPS,sizeof(unsigned) );
      
	unsigned INDEX,N,ID_J;
	for( unsigned i = 0 ; i < NHAPS ; ++i )
	  {
	    gzread( gzin,&INDEX,sizeof(unsigned) );
	    assert( rv.find(INDEX) == rv.end() );
	    gzread( gzin,&N,sizeof(unsigned) );
	    gamete_type g(N);
	    gzread( gzin,&N,sizeof(unsigned) );
	    for( unsigned j = 0 ; j < N ; ++j )
	      {
		gzread( gzin,&ID_J,sizeof(unsigned) );
		g.mutations.push_back(m[ID_J]);
	      }
	    gzread( gzin,&N,sizeof(unsigned) );
	    for( unsigned j = 0 ; j < N ; ++j )
	      {
		gzread( gzin,&ID_J,sizeof(unsigned) );
		g.smutations.push_back(m[ID_J]);
	      }
	    rv[INDEX] = gametes->insert(gametes->end(),g);
	  }
	return rv;
      }
    };
  }
}

#endif
