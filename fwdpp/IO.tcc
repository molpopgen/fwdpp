//  -*- C++ -*- 
#ifndef __FWDPP_IO_TCC__
#define __FWDPP_IO_TCC__

#include <fwdpp/forward_types.hpp>

#include <map>
#include <vector>
#include <algorithm>
#include <cassert>
#include <type_traits>

namespace KTfwd
{
  struct write_mutations
  {
    template< typename mutation_type,
	      typename list_type_allocator,
	      template<typename,typename> class list_type,
	      typename mutation_writer_type,
	      typename ostreamtype>
    std::pair< std::vector< typename list_type< mutation_type, list_type_allocator >::iterator >,
	       std::vector<unsigned> >
    operator()( list_type< mutation_type, list_type_allocator > * mutations,
		const mutation_writer_type & mw,
		ostreamtype & buffer) const
	      {
		typedef typename list_type< mutation_type, list_type_allocator >::iterator mlist_iterator;
		//initiate a list with the mutation information by iterator
		typedef std::vector< mlist_iterator > maptype;
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
      typename vector_type_allocator,
      template<typename,typename> class vector_type,
      typename ostreamtype>
    void operator()( const vector_type< gamete_type, vector_type_allocator > * gametes,
		     const std::vector< typename gamete_type::mutation_list_type_iterator > & mut_info,
		     const std::vector< unsigned > & indexes,
		     ostreamtype & buffer) const
      {
	typedef typename vector_type< gamete_type, vector_type_allocator>::const_iterator glist_iterator;
	unsigned NHAPS = unsigned(gametes->size());
	buffer.write( reinterpret_cast< char * >(&NHAPS), sizeof(unsigned) );
	for( glist_iterator gtr = gametes->begin() ; gtr != gametes->end() ; ++gtr )
	  {
	    unsigned N = gtr->n;
	    buffer.write( reinterpret_cast<char *>(&N), sizeof(unsigned) );
	    N = unsigned(gtr->mutations.size());
	    buffer.write( reinterpret_cast<char *>(&N), sizeof(unsigned) );
	    for( unsigned i = 0 ; i < N ; ++i )
	      {
		assert( std::find(mut_info.begin(),mut_info.end(),(gtr->mutations[i])) != mut_info.end() );
		unsigned INDEX = indexes[ std::vector<unsigned>::size_type(std::find(mut_info.begin(),mut_info.end(),(gtr->mutations[i])) - mut_info.begin()) ];
		buffer.write( reinterpret_cast< char * >(&INDEX), sizeof(unsigned) );
	      }
	    N = unsigned(gtr->smutations.size());
	    buffer.write( reinterpret_cast<char *>(&N), sizeof(unsigned) );
	    for( unsigned i = 0 ; i < N ; ++i )
	      {
		assert( std::find(mut_info.begin(),mut_info.end(),(gtr->smutations[i])) != mut_info.end() );
		unsigned INDEX = indexes[ std::vector<unsigned>::size_type(std::find(mut_info.begin(),mut_info.end(),(gtr->smutations[i])) - mut_info.begin()) ];
		buffer.write( reinterpret_cast< char * >(&INDEX), sizeof(unsigned) );
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
    std::map< unsigned, typename list_type<mutation_type,list_type_allocator>::iterator >
    operator()(list_type< mutation_type, list_type_allocator > * mutations,
	       const mutation_reader & mr,
	       istreamtype & in) const
    {
      typedef list_type< mutation_type, list_type_allocator > mlist;
      typedef std::map<unsigned,typename mlist::iterator> mut_info;
      
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
      typedef list_type< mutation_type, list_type_allocator > mlist;
      typedef std::map<unsigned,typename mlist::iterator> mut_info;
      
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
	      typename vector_type_allocator,
	      template<typename,typename> class vector_type,
	      typename istreamtype >
    void operator()(vector_type< gamete_type, vector_type_allocator > * gametes,
		    std::map<unsigned, typename gamete_type::mutation_list_type_iterator> & m,
		    istreamtype & in) const
    {
      unsigned NHAPS;
      in.read( reinterpret_cast< char * >(&NHAPS), sizeof(unsigned) );
      for( unsigned i = 0 ; i < NHAPS ; ++i)
	{
	  unsigned N;
	  in.read( reinterpret_cast< char * >(&N), sizeof(unsigned) );
	  gamete_type g(N);
	  in.read( reinterpret_cast< char * >(&N), sizeof(unsigned) );
	  for( unsigned j = 0 ; j < N ; ++j )
	    {
	      unsigned ID_J;
	      in.read( reinterpret_cast< char * >(&ID_J), sizeof(unsigned) );
	      g.mutations.push_back(m[ID_J]);
	    }
	  in.read( reinterpret_cast< char * >(&N), sizeof(unsigned) );
	  for( unsigned j = 0 ; j < N ; ++j )
	    {
	      unsigned ID_J;
	      in.read( reinterpret_cast< char * >(&ID_J), sizeof(unsigned) );
	      g.smutations.push_back(m[ID_J]);
	    }
	  gametes->push_back(g);
	}
    }

    template< typename gamete_type,
	      typename vector_type_allocator,
	      template<typename,typename> class vector_type>
    void operator()(vector_type< gamete_type, vector_type_allocator > * gametes,
		    std::map<unsigned, typename gamete_type::mutation_list_type_iterator> & m,
		    gzFile gzin) const
    {
      unsigned NHAPS;
      gzread( gzin, &NHAPS, sizeof(unsigned) );
      for( unsigned i = 0 ; i < NHAPS ; ++i)
	{
	  unsigned N;
	  gzread( gzin, &N, sizeof(unsigned) );
	  gamete_type g(N);
	  gzread( gzin, &N, sizeof(unsigned) );
	  for( unsigned j = 0 ; j < N ; ++j )
	    {
	      unsigned ID_J;
	      gzread( gzin, &ID_J, sizeof(unsigned) );
	      g.mutations.push_back(m[ID_J]);
	    }
	  gzread( gzin, &N, sizeof(unsigned) );
	  for( unsigned j = 0 ; j < N ; ++j )
	    {
	      unsigned ID_J;
	      gzread( gzin, &ID_J, sizeof(unsigned) );
	      g.smutations.push_back(m[ID_J]);
	    }
	  gametes->push_back(g);
	}
    }
  };

  template< typename gamete_type,
	    typename vector_type_allocator,
	    template<typename,typename> class vector_type,
	    typename mutation_type,
	    typename list_type_allocator,
	    template<typename,typename> class list_type,
	    typename mutation_writer_type,
	    typename ostreamtype>
  void write_binary_pop ( const vector_type< gamete_type, vector_type_allocator > * gametes,
			  list_type< mutation_type, list_type_allocator > * mutations,
			  const mutation_writer_type & mw,
			  ostreamtype & buffer)
  {	      
    static_assert( std::is_base_of<mutation_base,mutation_type>::value,
                   "mutation_type must inherit from KTfwd::mutation_base" );
    typedef gamete_base< typename gamete_type::mutation_type, typename gamete_type::mutation_list_type > gamete_base_type;
    static_assert( std::is_base_of<gamete_base_type,gamete_type>::value ||
                   std::is_same<gamete_base_type,gamete_type>::value,
                   "gamete_type must be, or inherit from, KTfwd::gamete_base<mutation_type,mutation_list_type>" );
    typedef typename list_type< mutation_type, list_type_allocator >::iterator mlist_iterator;
    //initiate a list with the mutation information by iterator
    typedef std::vector< mlist_iterator > maptype;
    std::pair< maptype, std::vector<unsigned> > mutdata = write_mutations()( mutations,mw,buffer );
    write_haplotypes()(gametes,mutdata.first,mutdata.second,buffer);
  }
  
  template< typename gamete_type,
	    typename vector_type_allocator,
	    template<typename,typename> class vector_type,
	    typename vector_type_allocator2,
	    template<typename,typename> class vector_type2,
	    typename mutation_type,
	    typename list_type_allocator,
	    template<typename,typename> class list_type,
	    typename mutation_writer_type,
	    typename ostreamtype>
  void write_binary_metapop ( const vector_type2<vector_type< gamete_type, vector_type_allocator >, vector_type_allocator2 > * gametes,
			      list_type< mutation_type, list_type_allocator > * mutations,
			      const mutation_writer_type & mw,
			      ostreamtype & buffer)
  {	      
    static_assert( std::is_base_of<mutation_base,mutation_type>::value,
                   "mutation_type must inherit from KTfwd::mutation_base" );
    typedef gamete_base< typename gamete_type::mutation_type, typename gamete_type::mutation_list_type > gamete_base_type;
    static_assert( std::is_base_of<gamete_base_type,gamete_type>::value ||
                   std::is_same<gamete_base_type,gamete_type>::value,
                   "gamete_type must be, or inherit from, KTfwd::gamete_base<mutation_type,mutation_list_type>" );
    typedef typename list_type< mutation_type, list_type_allocator >::iterator mlist_iterator;
    
    //initiate a list with the mutation information by iterator
    typedef std::vector< mlist_iterator > maptype;
    std::pair< maptype, std::vector<unsigned> > mutdata = write_mutations()( mutations,mw,buffer );
    unsigned NPOPS = gametes->size();
    buffer.write( reinterpret_cast< char * >(&NPOPS), sizeof(unsigned) );
    for( typename  vector_type2<vector_type< gamete_type, vector_type_allocator >, 
	   vector_type_allocator2>::const_iterator pop = gametes->begin() ; 
	 pop < gametes->end() ; ++pop )
      {
	write_haplotypes()(&*pop,mutdata.first,mutdata.second,buffer);
      }
  }

  template< typename gamete_type,
	    typename vector_type_allocator,
	    template<typename,typename> class vector_type,
	    typename mutation_type,
	    typename list_type_allocator,
	    template<typename,typename> class list_type,
	    typename mutation_reader,
	    typename istreamtype>
  void read_binary_pop ( vector_type< gamete_type, vector_type_allocator > * gametes,
			 list_type< mutation_type, list_type_allocator > * mutations,
			 const mutation_reader & mr,
			 istreamtype & in)
  {
    gametes->clear();
    mutations->clear();
    
    static_assert( std::is_base_of<mutation_base,mutation_type>::value,
                   "mutation_type must inherit from KTfwd::mutation_base" );
    typedef gamete_base< typename gamete_type::mutation_type, typename gamete_type::mutation_list_type > gamete_base_type;
    static_assert( std::is_base_of<gamete_base_type,gamete_type>::value ||
                   std::is_same<gamete_base_type,gamete_type>::value,
                   "gamete_type must be, or inherit from, KTfwd::gamete_base<mutation_type,mutation_list_type>" );
    typedef list_type< mutation_type, list_type_allocator > mlist;
    typedef std::map<unsigned,typename mlist::iterator> mut_info;
    
    mut_info m = read_mutations()(mutations,mr,in);
    read_haplotypes()(gametes,m,in);
  }
			 
  template< typename gamete_type,
	    typename vector_type_allocator,
	    template<typename,typename> class vector_type,
	    typename vector_type_allocator2,
	    template<typename,typename> class vector_type2,
	    typename mutation_type,
	    typename list_type_allocator,
	    template<typename,typename> class list_type,
	    typename mutation_reader,
	    typename istreamtype>
  void read_binary_metapop ( vector_type2< vector_type< gamete_type, vector_type_allocator >, vector_type_allocator2 > * gametes,
			     list_type< mutation_type, list_type_allocator > * mutations,
			     const mutation_reader & mr,
			     istreamtype & in)
  {
    gametes->clear();
    mutations->clear();
    
    static_assert( std::is_base_of<mutation_base,mutation_type>::value,
                   "mutation_type must inherit from KTfwd::mutation_base" );
    typedef gamete_base< typename gamete_type::mutation_type, typename gamete_type::mutation_list_type > gamete_base_type;
    static_assert( std::is_base_of<gamete_base_type,gamete_type>::value ||
                   std::is_same<gamete_base_type,gamete_type>::value,
                   "gamete_type must be, or inherit from, KTfwd::gamete_base<mutation_type,mutation_list_type>" );
    typedef list_type< mutation_type, list_type_allocator > mlist;
    typedef std::map<unsigned,typename mlist::iterator> mut_info;
    
    mut_info m = read_mutations()(mutations,mr,in);
    unsigned NPOP = 0;
    in.read( reinterpret_cast< char * >(&NPOP), sizeof(unsigned) );
    gametes->resize(NPOP);
    for(typename vector_type2< vector_type< gamete_type, vector_type_allocator >, 
	  vector_type_allocator2 >::iterator pop = gametes->begin() ; 
	pop < gametes->end() ; ++pop)
      {
	pop->clear();
	read_haplotypes()(&*pop,m,in);
      }
  }

  template< typename gamete_type,
	    typename vector_type_allocator,
	    template<typename,typename> class vector_type,
	    typename vector_type_allocator2,
	    template<typename,typename> class vector_type2,
	    typename mutation_type,
	    typename list_type_allocator,
	    template<typename,typename> class list_type,
	    typename mutation_reader>
  void read_binary_metapop ( vector_type2< vector_type< gamete_type, vector_type_allocator >, vector_type_allocator2 > * gametes,
			     list_type< mutation_type, list_type_allocator > * mutations,
			     const mutation_reader & mr,
			     gzFile gzin )
  {
    gametes->clear();
    mutations->clear();
    static_assert( std::is_base_of<mutation_base,mutation_type>::value,
                   "mutation_type must inherit from KTfwd::mutation_base" );    
    typedef gamete_base< typename gamete_type::mutation_type, typename gamete_type::mutation_list_type > gamete_base_type;
    static_assert( std::is_base_of<gamete_base_type,gamete_type>::value ||
                   std::is_same<gamete_base_type,gamete_type>::value,
                   "gamete_type must be, or inherit from, KTfwd::gamete_base<mutation_type,mutation_list_type>" );
    typedef list_type< mutation_type, list_type_allocator > mlist;
    typedef std::map<unsigned,typename mlist::iterator> mut_info;
    
    mut_info m = read_mutations()(mutations,mr,gzin);
    unsigned NPOP = 0;
    gzread(gzin,&NPOP,sizeof(unsigned));
    gametes->resize(NPOP);
    for(typename vector_type2< vector_type< gamete_type, vector_type_allocator >, 
	  vector_type_allocator2 >::iterator pop = gametes->begin() ; 
	pop < gametes->end() ; ++pop)
      {
	pop->clear();
	read_haplotypes()(&*pop,m,gzin);
      }
  }
  
  //Binary I/O for individual-based simulation

  struct write_haplotypes_ind
  {
    template<
      typename gamete_type,
      typename list_type_allocator,
      template<typename,typename> class list_type,
      typename ostreamtype>
    std::pair< std::vector< typename list_type< gamete_type, list_type_allocator >::const_iterator >,
	       std::vector<unsigned> >
    operator()( const list_type< gamete_type, list_type_allocator > * gametes,
		const std::vector< typename gamete_type::mutation_list_type_iterator > & mut_info,
		const std::vector< unsigned > & indexes,
		ostreamtype & buffer) const
    {
      typedef typename list_type< gamete_type, list_type_allocator>::const_iterator glist_iterator;
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

  template< typename gamete_type,
	    typename gamete_list_type_allocator,
	    template<typename,typename> class gamete_list_type,
	    typename mutation_type,
	    typename mutation_list_type_allocator,
	    template<typename,typename> class mutation_list_type,
	    typename vector_type_allocator,
	    template<typename,typename> class diploid_vector_type,
	    typename mutation_writer_type,
	    typename ostreamtype>
  void write_binary_pop ( const gamete_list_type< gamete_type, gamete_list_type_allocator > * gametes,
			  mutation_list_type< mutation_type, mutation_list_type_allocator > * mutations,
			  const diploid_vector_type< std::pair< typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator,
								typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator >,
						     vector_type_allocator > * diploids,
			  const mutation_writer_type & mw,
			  ostreamtype & buffer)
  {
    typedef typename mutation_list_type< mutation_type, mutation_list_type_allocator >::iterator mlist_iterator;
    typedef typename gamete_list_type< gamete_type, gamete_list_type_allocator >::const_iterator glist_iterator;
    //initiate a list with the mutation information by iterator
    typedef std::vector< mlist_iterator > maptype;
    std::pair< maptype, std::vector<unsigned> > mutdata = write_mutations()( mutations,mw,buffer ); 
    
    typedef std::vector< glist_iterator > gmaptype;
    
    std::pair<gmaptype, std::vector<unsigned> > gamdata = write_haplotypes_ind()( gametes, mutdata.first, mutdata.second, buffer );
    
    //the diploids are now trivial to write
    typedef typename  diploid_vector_type< std::pair< typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator,
						      typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator >,
					   vector_type_allocator >::const_iterator dptr;

    unsigned NDIPS = unsigned(diploids->size());
    buffer.write( reinterpret_cast<char *>(&NDIPS), sizeof(unsigned) );
    for( dptr dip = diploids->begin() ; dip != diploids->end() ; ++dip )
      {
	unsigned c = gamdata.second[ std::vector<unsigned>::size_type(std::find( gamdata.first.begin(),gamdata.first.end(),dip->first ) - gamdata.first.begin()) ];
	buffer.write( reinterpret_cast<char *>(&c), sizeof(unsigned) );
	c = gamdata.second[ std::vector<unsigned>::size_type(std::find( gamdata.first.begin(),gamdata.first.end(),dip->second ) - gamdata.first.begin()) ];
	buffer.write( reinterpret_cast<char *>(&c), sizeof(unsigned) );
      }
  }

  struct read_haplotypes_ind
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

  template< typename gamete_type,
	    typename gamete_list_type_allocator,
	    template<typename,typename> class gamete_list_type,
	    typename mutation_type,
	    typename mutation_list_type_allocator,
	    template<typename,typename> class mutation_list_type,
	    typename vector_type_allocator,
	    template<typename,typename> class diploid_vector_type,
	    typename mutation_reader_type,
	    typename istreamtype>
  void read_binary_pop (  gamete_list_type< gamete_type, gamete_list_type_allocator > * gametes,
			  mutation_list_type< mutation_type, mutation_list_type_allocator > * mutations,
			  diploid_vector_type< std::pair< typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator,
							  typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator >,
					       vector_type_allocator > * diploids,
			  const mutation_reader_type & mr,
			  istreamtype & in)
  {
    gametes->clear();
    mutations->clear();
    diploids->clear();
    
    typedef mutation_list_type< mutation_type, mutation_list_type_allocator > mlist;
    typedef std::map<unsigned,typename mlist::iterator> mut_info;
    
    mut_info m = read_mutations()(mutations,mr,in); 
    
    typedef gamete_list_type< gamete_type, gamete_list_type_allocator > glist;
    typedef std::map<unsigned,typename glist::iterator> gam_info;
    
    gam_info g = read_haplotypes_ind()(gametes,m,in);
    
    unsigned NDIPS,c;
    in.read( reinterpret_cast<char *>(&NDIPS), sizeof(unsigned) );
    
    typedef typename  diploid_vector_type< std::pair< typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator,
						      typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator >,
					   vector_type_allocator >::iterator dptr;

    diploids->resize(NDIPS);
    for( dptr dp = diploids->begin() ; dp != diploids->end() ; ++dp )
      {
	in.read( reinterpret_cast<char *>(&c), sizeof(unsigned) );
	dp->first = g[c];
	in.read( reinterpret_cast<char *>(&c), sizeof(unsigned) );
	dp->second = g[c];
      }
  }
  
  template< typename gamete_type,
	    typename gamete_list_type_allocator,
	    template<typename,typename> class gamete_list_type,
	    typename mutation_type,
	    typename mutation_list_type_allocator,
	    template<typename,typename> class mutation_list_type,
	    typename vector_type_allocator,
	    template<typename,typename> class diploid_vector_type,
	    typename mutation_reader_type>
  void read_binary_pop (  gamete_list_type< gamete_type, gamete_list_type_allocator > * gametes,
			  mutation_list_type< mutation_type, mutation_list_type_allocator > * mutations,
			  diploid_vector_type< std::pair< typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator,
							  typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator >,
					       vector_type_allocator > * diploids,
			  const mutation_reader_type & mr,
			  gzFile gzin )
  {
    gametes->clear();
    mutations->clear();
    diploids->clear();
    
    typedef mutation_list_type< mutation_type, mutation_list_type_allocator > mlist;
    typedef std::map<unsigned,typename mlist::iterator> mut_info;
    
    mut_info m = read_mutations()(mutations,mr,gzin); 
    
    typedef gamete_list_type< gamete_type, gamete_list_type_allocator > glist;
    typedef std::map<unsigned,typename glist::iterator> gam_info;
    
    gam_info g = read_haplotypes_ind()(gametes,m,gzin);
    
    unsigned NDIPS,c;
    gzread( gzin,&NDIPS,sizeof(unsigned) );
    
    typedef typename  diploid_vector_type< std::pair< typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator,
						      typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator >,
					   vector_type_allocator >::iterator dptr;

    diploids->resize(NDIPS);
    for( dptr dp = diploids->begin() ; dp != diploids->end() ; ++dp )
      {
	gzread( gzin,&c,sizeof(unsigned) );
	dp->first = g[c];
	gzread( gzin,&c,sizeof(unsigned) );
	dp->second = g[c];
      }
  }


  template< typename gamete_type,
	    typename gamete_list_type_allocator,
	    template<typename,typename> class gamete_list_type,
	    typename metapop_vector_type_allocator,
	    template<typename,typename> class metapop_vector_type,
	    typename mutation_type,
	    typename mutation_list_type_allocator,
	    template<typename,typename> class mutation_list_type,
	    typename vector_type_allocator,
	    template<typename,typename> class diploid_vector_type,
	    typename diploid_vv_type_allocator,
	    template<typename,typename> class diploid_vv_type,
	    typename mutation_writer_type,
	    typename ostreamtype>
  void write_binary_metapop ( const metapop_vector_type< gamete_list_type< gamete_type, gamete_list_type_allocator >,
			      metapop_vector_type_allocator> * metapop,
			      mutation_list_type< mutation_type, mutation_list_type_allocator > * mutations,
			      const diploid_vv_type < diploid_vector_type< std::pair< typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator,
										      typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator >,
									   vector_type_allocator >,
			      diploid_vv_type_allocator > * diploids,
			      const mutation_writer_type & mw,
			      ostreamtype & buffer)
  {
    unsigned NPOP = unsigned(metapop->size());
    
    buffer.write( reinterpret_cast<char *>(&NPOP), sizeof(unsigned) );
    
    typedef mutation_list_type< mutation_type, mutation_list_type_allocator > mlist;
    typedef typename mlist::iterator mlist_iterator;
    
    typedef std::vector< mlist_iterator > maptype;
    std::pair< maptype, std::vector<unsigned> > mutdata = write_mutations()( mutations,mw,buffer ); 
    
    typename metapop_vector_type< gamete_list_type< gamete_type, gamete_list_type_allocator >,
				  metapop_vector_type_allocator>::const_iterator pop_ptr = metapop->begin();
  
  typename diploid_vv_type < diploid_vector_type< std::pair< typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator,
							     typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator >,
						  vector_type_allocator >,
			     diploid_vv_type_allocator >::const_iterator dip_ptr = diploids->begin();
  
  typedef typename diploid_vector_type< std::pair< typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator,
						   typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator >,
					vector_type_allocator >::const_iterator dptr;
  
  typedef std::vector< typename gamete_list_type< gamete_type, gamete_list_type_allocator >::const_iterator > gmaptype;
  
  unsigned NDIPS,c;
  for( unsigned pop = 0 ; pop < NPOP ; ++pop,++pop_ptr,++dip_ptr )
    {
      std::pair<gmaptype, std::vector<unsigned> > gamdata = write_haplotypes_ind()( &*pop_ptr, mutdata.first, mutdata.second, buffer );
      
      NDIPS = unsigned(dip_ptr->size());
      
      buffer.write( reinterpret_cast<char *>(&NDIPS),sizeof(unsigned) );
      
      for( dptr dip = dip_ptr->begin() ; dip != dip_ptr->end() ; ++dip )
	{
	  c = gamdata.second[ std::vector<unsigned>::size_type(std::find( gamdata.first.begin(), gamdata.first.end(), dip->first ) - gamdata.first.begin()) ];
	  buffer.write( reinterpret_cast<char *>(&c),sizeof(unsigned) );
	  c = gamdata.second[ std::vector<unsigned>::size_type(std::find( gamdata.first.begin(), gamdata.first.end(), dip->second ) - gamdata.first.begin()) ];
	  buffer.write( reinterpret_cast<char *>(&c),sizeof(unsigned) );
	}
    }
  }
  
  template< typename gamete_type,
	    typename gamete_list_type_allocator,
	    template<typename,typename> class gamete_list_type,
	    typename metapop_vector_type_allocator,
	    template<typename,typename> class metapop_vector_type,
	    typename mutation_type,
	    typename mutation_list_type_allocator,
	    template<typename,typename> class mutation_list_type,
	    typename vector_type_allocator,
	    template<typename,typename> class diploid_vector_type,
	    typename diploid_vv_type_allocator,
	    template<typename,typename> class diploid_vv_type,
	    typename mutation_reader_type,
	    typename istreamtype>
  void read_binary_metapop ( metapop_vector_type< gamete_list_type< gamete_type, gamete_list_type_allocator >,
			     metapop_vector_type_allocator> * metapop,
			     mutation_list_type< mutation_type, mutation_list_type_allocator > * mutations,
			     diploid_vv_type < diploid_vector_type< std::pair< typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator,
									       typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator >,
								    vector_type_allocator >,
			     diploid_vv_type_allocator > * diploids,
			     const mutation_reader_type & mr,
			     istreamtype & in)
  {
    metapop->clear();
    mutations->clear();
    diploids->clear();
    
    typedef mutation_list_type< mutation_type, mutation_list_type_allocator > mlist;
    typedef std::map<unsigned,typename mlist::iterator> mut_info;
    
    unsigned NPOP;
    in.read( reinterpret_cast< char * >(&NPOP), sizeof(unsigned) );
    
    mut_info m = read_mutations()(mutations,mr,in); 
    
    metapop->resize(NPOP);
    diploids->resize(NPOP);

    typename metapop_vector_type< gamete_list_type< gamete_type, gamete_list_type_allocator >,
				  metapop_vector_type_allocator>::iterator pop_ptr = metapop->begin();
  
    typename diploid_vv_type < diploid_vector_type< std::pair< typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator,
							       typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator >,
						    vector_type_allocator >,
			       diploid_vv_type_allocator >::iterator dip_ptr = diploids->begin();
    
    typedef typename diploid_vector_type< std::pair< typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator,
						     typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator >,
					  vector_type_allocator >::iterator dptr;

typedef gamete_list_type< gamete_type, gamete_list_type_allocator > glist;
typedef std::map<unsigned,typename glist::iterator> gam_info;
  
  unsigned NDIPS,c;
  for( unsigned pop=0 ; pop < NPOP ; ++pop,++pop_ptr,++dip_ptr )
    {
      
      gam_info g = read_haplotypes_ind()(&*pop_ptr,m,in);	
      
      in.read(reinterpret_cast<char*>(&NDIPS),sizeof(unsigned));
      
      dip_ptr->resize(NDIPS);
      
      for( dptr dip = dip_ptr->begin() ; dip != dip_ptr->end() ; ++dip )
	{
	  in.read(reinterpret_cast<char *>(&c),sizeof(unsigned));
	  dip->first = g[c];
	  in.read(reinterpret_cast<char *>(&c),sizeof(unsigned));
	  dip->second = g[c];
	}
    }
  }

  template< typename gamete_type,
	    typename gamete_list_type_allocator,
	    template<typename,typename> class gamete_list_type,
	    typename metapop_vector_type_allocator,
	    template<typename,typename> class metapop_vector_type,
	    typename mutation_type,
	    typename mutation_list_type_allocator,
	    template<typename,typename> class mutation_list_type,
	    typename vector_type_allocator,
	    template<typename,typename> class diploid_vector_type,
	    typename diploid_vv_type_allocator,
	    template<typename,typename> class diploid_vv_type,
	    typename mutation_reader_type>
  void read_binary_metapop ( metapop_vector_type< gamete_list_type< gamete_type, gamete_list_type_allocator >,
			     metapop_vector_type_allocator> * metapop,
			     mutation_list_type< mutation_type, mutation_list_type_allocator > * mutations,
			     diploid_vv_type < diploid_vector_type< std::pair< typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator,
									       typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator >,
								    vector_type_allocator >,
			     diploid_vv_type_allocator > * diploids,
			     const mutation_reader_type & mr,
			     gzFile gzin)
  {
    metapop->clear();
    mutations->clear();
    diploids->clear();
    
    typedef mutation_list_type< mutation_type, mutation_list_type_allocator > mlist;
    typedef std::map<unsigned,typename mlist::iterator> mut_info;
    
    unsigned NPOP;
    gzread( gzin,&NPOP,sizeof(unsigned) );
    
    mut_info m = read_mutations()(mutations,mr,gzin); 
    
    metapop->resize(NPOP);
    diploids->resize(NPOP);

    typename metapop_vector_type< gamete_list_type< gamete_type, gamete_list_type_allocator >,
				  metapop_vector_type_allocator>::iterator pop_ptr = metapop->begin();
  
    typename diploid_vv_type < diploid_vector_type< std::pair< typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator,
							       typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator >,
						    vector_type_allocator >,
			       diploid_vv_type_allocator >::iterator dip_ptr = diploids->begin();
    
    typedef typename diploid_vector_type< std::pair< typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator,
						     typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator >,
					  vector_type_allocator >::iterator dptr;

typedef gamete_list_type< gamete_type, gamete_list_type_allocator > glist;
typedef std::map<unsigned,typename glist::iterator> gam_info;
  
  unsigned NDIPS,c;
  for( unsigned pop=0 ; pop < NPOP ; ++pop,++pop_ptr,++dip_ptr )
    {
      
      gam_info g = read_haplotypes_ind()(&*pop_ptr,m,gzin);	
      
      gzread( gzin,&NDIPS,sizeof(unsigned) );
      dip_ptr->resize(NDIPS);
      
      for( dptr dip = dip_ptr->begin() ; dip != dip_ptr->end() ; ++dip )
	{
	  gzread( gzin,&c,sizeof(unsigned) );
	  dip->first = g[c];
	  gzread( gzin,&c,sizeof(unsigned) );
	  dip->second = g[c];
	}
    }
  }
}//ns KTfwd

#endif
