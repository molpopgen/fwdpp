#ifndef __FWDPP_SUGAR_SINGLEPOP_SERIALIZATON_HPP__
#define __FWDPP_SUGAR_SINGLEPOP_SERIALIZATON_HPP__

#include <iosfwd>
#include <sstream>
#include <algorithm>
#include <type_traits>
#include <utility>
#include <functional>
#include <fwdpp/IO.hpp>
#include <fwdpp/sugar/poptypes/tags.hpp>
#include <fwdpp/sugar/popgenmut.hpp>
#include <fwdpp/sugar/generalmut.hpp>

namespace KTfwd
{
  /*!
    \brief Facilitates serialization of mutation
    types supported by the fwdpp sugar library
    \ingroup sugar
  */
  struct mutation_writer
  {
    /*!
      \brief overload for KTfwd::popgenmut and ostreams
     */
    using result_type = void;
    template<typename mutation_t>
    inline typename std::enable_if<std::is_same<mutation_t,popgenmut>::value,result_type>::type
    operator()( const mutation_t &m,
		std::ostream & buffer) const
    {
      buffer.write( reinterpret_cast<const char *>(&m.n),sizeof(unsigned));
      buffer.write( reinterpret_cast<const char *>(&m.g),sizeof(unsigned));
      buffer.write( reinterpret_cast<const char *>(&m.pos),sizeof(double));
      buffer.write( reinterpret_cast<const char *>(&m.s),sizeof(double));
      buffer.write( reinterpret_cast<const char *>(&m.h),sizeof(double));
    }

    /*!
      \brief overload for KTfwd::popgenmut and zlib/gzFile
    */
    template<typename mutation_t>
    inline typename std::enable_if<std::is_same<mutation_t,popgenmut>::value,result_type>::type
    operator()( const mutation_t &m,
		gzFile gzout) const
    {
      gzwrite(gzout, reinterpret_cast<const char *>(&m.n),sizeof(unsigned));
      gzwrite(gzout, reinterpret_cast<const char *>(&m.g),sizeof(unsigned));
      gzwrite(gzout, reinterpret_cast<const char *>(&m.pos),sizeof(double));
      gzwrite(gzout, reinterpret_cast<const char *>(&m.s),sizeof(double));
      gzwrite(gzout, reinterpret_cast<const char *>(&m.h),sizeof(double));
    }
    
    /*!
      \brief overload for KTfwd::mutation and ostreams
     */
    template<typename mutation_t>
    inline typename std::enable_if<std::is_same<mutation_t,mutation>::value,result_type>::type
    operator()( const mutation_t &m,
		std::ostream & buffer) const
    {
      buffer.write( reinterpret_cast<const char *>(&m.n),sizeof(unsigned));
      buffer.write( reinterpret_cast<const char *>(&m.pos),sizeof(double));
      buffer.write( reinterpret_cast<const char *>(&m.s),sizeof(double));
      buffer.write( reinterpret_cast<const char *>(&m.h),sizeof(double));
    }

    /*!
      \brief overload for KTfwd::mutation and zlib/gzFile
    */
    template<typename mutation_t>
    inline typename std::enable_if<std::is_same<mutation_t,mutation>::value,result_type>::type
    operator()( const mutation_t &m,
		gzFile gzout) const
    {
      gzwrite(gzout, reinterpret_cast<const char *>(&m.n),sizeof(unsigned));
      gzwrite(gzout, reinterpret_cast<const char *>(&m.pos),sizeof(double));
      gzwrite(gzout, reinterpret_cast<const char *>(&m.s),sizeof(double));
      gzwrite(gzout, reinterpret_cast<const char *>(&m.h),sizeof(double));
    }

    //! \brief overload for KTfwd::generalmut and ostream
    template<typename mutation_t,
	     std::size_t N = std::tuple_size<typename mutation_t::array_t>()>
    inline typename std::enable_if<std::is_same<mutation_t,generalmut<N> >::value,result_type>::type
    operator()(const mutation_t & t, std::ostream & buffer)
    {
      buffer.write( reinterpret_cast<const char *>(&t.n),sizeof(unsigned));
      buffer.write( reinterpret_cast<const char *>(&t.g),sizeof(unsigned));
      buffer.write( reinterpret_cast<const char *>(&t.pos),sizeof(double));
      //Write mutation types
      buffer.write( reinterpret_cast<const char *>(&t.s[0]),N*sizeof(double));
      buffer.write( reinterpret_cast<const char *>(&t.h[0]),N*sizeof(double));
    }

    //! \brief overload for KTfwd::generalmut and zlib/gzFile
    template<typename mutation_t,
	     std::size_t N = std::tuple_size<typename mutation_t::array_t>()>
    inline typename std::enable_if<std::is_same<mutation_t,generalmut<N> >::value,result_type>::type
    operator()(const mutation_t & t, gzFile gzout)
    {
      gzwrite(gzout, reinterpret_cast<const char *>(&t.n),sizeof(unsigned));
      gzwrite(gzout, reinterpret_cast<const char *>(&t.g),sizeof(unsigned));
      gzwrite(gzout, reinterpret_cast<const char *>(&t.pos),sizeof(double));
      //Write mutation types
      gzwrite(gzout, reinterpret_cast<const char *>(&t.s[0]),N*sizeof(double));
      gzwrite(gzout, reinterpret_cast<const char *>(&t.h[0]),N*sizeof(double));
    }
  };

  /*!
    \brief Facilitates serialization of mutation
    types supported by the fwdpp sugar library.

    The template parameter must be derived from
    KTfwd::mutation_base

    \ingroup sugar
  */
  template<typename mutation_t>
  struct mutation_reader
  {
    //! The return value of operator()
    using result_type = mutation_t;
    /*!
      \brief overload for KTfwd::popgenmut and istreams
     */
    template<typename U = mutation_t>
    inline typename std::enable_if<std::is_same<U,popgenmut>::value,result_type>::type
    operator()( std::istream & in ) const
    {
      unsigned n,g;
      double pos,s,h;
      in.read( reinterpret_cast<char *>(&n),sizeof(unsigned));
      in.read( reinterpret_cast<char *>(&g),sizeof(unsigned));
      in.read( reinterpret_cast<char *>(&pos),sizeof(double));
      in.read( reinterpret_cast<char *>(&s),sizeof(double));
      in.read( reinterpret_cast<char *>(&h),sizeof(double));
      return result_type(pos,s,h,g,n);
    }

    /*!
      \brief overload for KTfwd::popgenmut and zlib/gzFile
    */
    template<typename U = mutation_t>
    inline typename std::enable_if<std::is_same<U,popgenmut>::value,result_type>::type
    operator()( gzFile in ) const
    {
      unsigned n,g;
      double pos,s,h;
      gzread(in,&n,sizeof(unsigned));
      gzread(in,&g,sizeof(unsigned));
      gzread(in,&pos,sizeof(double));
      gzread(in,&s,sizeof(double));
      gzread(in,&h,sizeof(double));
      return result_type(pos,s,h,g,n);
    }
    
    /*!
      \brief overload for KTfwd::mutation and istreams
     */
    template<typename U = mutation_t>
    inline typename std::enable_if<std::is_same<U,mutation>::value,result_type>::type
    operator()( std::istream & in ) const
    {
      unsigned n;
      double pos,s,h;
      in.read( reinterpret_cast<char *>(&n),sizeof(unsigned));
      in.read( reinterpret_cast<char *>(&pos),sizeof(double));
      in.read( reinterpret_cast<char *>(&s),sizeof(double));
      in.read( reinterpret_cast<char *>(&h),sizeof(double));
      return result_type(pos,s,n,h);
    }

        
    /*!
      \brief overload for KTfwd::mutation and gzFile
    */
    template<typename U = mutation_t>
    inline typename std::enable_if<std::is_same<U,mutation>::value,result_type>::type
    operator()( gzFile in ) const
    {
      unsigned n;
      double pos,s,h;
      gzread(in,&n,sizeof(unsigned));
      gzread(in,&pos,sizeof(double));
      gzread(in,&s,sizeof(double));
      gzread(in,&h,sizeof(double));
      return result_type(pos,s,n,h);
    }

    //! \brief overalod for KTfwd::genericmut and std::istream
    template<typename U = mutation_t>
    inline typename std::enable_if<std::is_same<mutation_t,generalmut<std::tuple_size<typename U::array_t>::value> >::value,result_type>::type
    operator()(std::istream & buffer)
    {
      unsigned n,g;
      double pos;
      std::array<double,std::tuple_size<typename U::array_t>::value> s,h;
      buffer.read( reinterpret_cast<char *>(&n),sizeof(unsigned));
      buffer.read( reinterpret_cast<char *>(&g),sizeof(unsigned));
      buffer.read( reinterpret_cast<char *>(&pos),sizeof(double));
      //Write mutation types
      buffer.read( reinterpret_cast<char *>(&s[0]),std::tuple_size<typename U::array_t>::value*sizeof(double));
      buffer.read( reinterpret_cast<char *>(&h[0]),std::tuple_size<typename U::array_t>::value*sizeof(double));
      return generalmut<std::tuple_size<typename U::array_t>::value>(s,h,pos,n,g);
    }

    //! \brief overalod for KTfwd::genericmut and zlib/gzFile
    template<typename U = mutation_t>
    inline typename std::enable_if<std::is_same<mutation_t,generalmut<std::tuple_size<typename U::array_t>::value> >::value,result_type>::type
    operator()(gzFile in)
    {
      unsigned n,g;
      double pos;
      std::array<double,std::tuple_size<typename U::array_t>::value> s,h;
      gzread(in, &n,sizeof(unsigned));
      gzread(in, &g,sizeof(unsigned));
      gzread(in, &pos,sizeof(double));
      //Write mutation types
      gzread(in,&s[0],std::tuple_size<typename U::array_t>::value*sizeof(double));
      gzread(in,&h[0],std::tuple_size<typename U::array_t>::value*sizeof(double));
      return generalmut<std::tuple_size<typename U::array_t>::value>(s,h,pos,n,g);
    }
  };

  /*!
    \brief Serialize populations.
    \ingroup sugar
   */
  struct serialize
  {
    using result_type = void;
    mutable std::stringstream buffer;

    //!Default constructor
    serialize() : buffer(std::string()) {}

    //! Move constructor.  Req'd for this to be a member type of another class
    serialize( serialize && __s) : buffer(__s.buffer.str())
    {
    }
    
    /*!
      \brief Overload for single population simulations
     */
    template<typename sugarpop_t,
	     typename writer_t,
	     typename diploid_writer_t = diploidIOplaceholder>
    inline typename std::enable_if<std::is_same<typename sugarpop_t::popmodel_t,sugar::SINGLEPOP_TAG>::value,result_type>::type
    operator()( const sugarpop_t & pop,
		const writer_t & wt,
		const diploid_writer_t & dw = diploid_writer_t()) const
    {
      buffer.write( reinterpret_cast<const char*>(&pop.N), sizeof(unsigned) );
      write_binary_pop(&pop.gametes,&pop.mutations,&pop.diploids,wt,buffer,dw);
      //Step 2: output fixations 
      unsigned temp = pop.fixations.size();
      buffer.write( reinterpret_cast<char*>(&temp), sizeof(unsigned) );
      if( temp )
	{
	  std::for_each( pop.fixations.begin(), pop.fixations.end(),
			 std::bind(std::cref(wt),std::placeholders::_1,std::ref(buffer)) );
	  //Step 3:the fixation times
	  buffer.write( reinterpret_cast<const char *>(&pop.fixation_times[0]), pop.fixation_times.size()*sizeof(unsigned) );
	}
    }

    template<typename sugarpop_t,
	     typename writer_t,
	     typename diploid_writer_t = diploidIOplaceholder>
    inline typename std::enable_if<std::is_same<typename sugarpop_t::popmodel_t,sugar::MULTILOCPOP_TAG>::value,result_type>::type
    operator()( const sugarpop_t & pop,
		const writer_t & wt,
		const diploid_writer_t & dw = diploid_writer_t()) const
    {
      buffer.write( reinterpret_cast<const char*>(&pop.N), sizeof(unsigned) );
      write_binary_pop_mloc(&pop.gametes,&pop.mutations,&pop.diploids,wt,buffer,dw);
      //Step 2: output fixations 
      unsigned temp = pop.fixations.size();
      buffer.write( reinterpret_cast<char*>(&temp), sizeof(unsigned) );
      if( temp )
	{
	  std::for_each( pop.fixations.begin(), pop.fixations.end(),
			 std::bind(std::cref(wt),std::placeholders::_1,std::ref(buffer)) );
	  //Step 3:the fixation times
	  buffer.write( reinterpret_cast<const char *>(&pop.fixation_times[0]), pop.fixation_times.size()*sizeof(unsigned) );
	}
    }

    /*!
      \brief Overload for metapopulation simulations
     */
    template<typename sugarpop_t,
	     typename writer_t,
	     typename diploid_writer_t = diploidIOplaceholder>
    inline typename std::enable_if<std::is_same<typename sugarpop_t::popmodel_t,sugar::METAPOP_TAG>::value,result_type>::type
    operator()( const sugarpop_t & pop,
		const writer_t & wt,
		const diploid_writer_t & dw = diploid_writer_t()) const
    {
      unsigned npops = pop.Ns.size();
      buffer.write(reinterpret_cast<char *>(&npops),sizeof(unsigned));
      buffer.write( reinterpret_cast<const char*>(&pop.Ns[0]), npops*sizeof(unsigned) );
      write_binary_metapop(&pop.gametes,&pop.mutations,&pop.diploids,wt,buffer,dw);
      //Step 2: output fixations 
      unsigned temp = pop.fixations.size();
      buffer.write( reinterpret_cast<char*>(&temp), sizeof(unsigned) );
      if( temp )
	{
	  std::for_each( pop.fixations.begin(), pop.fixations.end(),
			 std::bind(std::cref(wt),std::placeholders::_1,std::ref(buffer)) );
	  //Step 3:the fixation times
	  buffer.write( reinterpret_cast<const char *>(&pop.fixation_times[0]), pop.fixation_times.size()*sizeof(unsigned) );
	}
    }
  };


  /*!
    \brief Deserialize population objects
    \ingroup sugar
   */
  struct deserialize
  {
    //! The return type for operator()
    using result_type = void;
    /*!
      \brief Overload for single population simulations
     */
    template<typename sugarpop_t,
	     typename reader_t,
	     typename diploid_reader_t = diploidIOplaceholder>
    inline typename std::enable_if<std::is_same<typename sugarpop_t::popmodel_t,sugar::SINGLEPOP_TAG>::value,result_type>::type
    operator()( sugarpop_t & pop,
		const serialize & s,
		const reader_t & rt,
		const diploid_reader_t & dr = diploid_reader_t()) const
    {
      pop.clear();
      //Step 0: read N
      s.buffer.read( reinterpret_cast<char*>(&pop.N),sizeof(unsigned) );
      KTfwd::read_binary_pop( &pop.gametes,&pop.mutations,&pop.diploids,rt,s.buffer,dr );
      unsigned temp;
      s.buffer.read( reinterpret_cast<char*>(&temp),sizeof(unsigned) );
      for( unsigned m=0;m<temp ;++m )
    	{
    	  typename reader_t::result_type mm = rt(s.buffer);
    	  pop.fixations.emplace_back( std::move(mm) );
    	}
      pop.fixation_times.resize(temp);
      if(temp)
	{
	  s.buffer.read( reinterpret_cast<char*>(&pop.fixation_times[0]), temp*sizeof(unsigned) );
	}
      s.buffer.seekg(0);
      //Finally, fill the lookup table:
      std::for_each( pop.mutations.begin(), pop.mutations.end(),
    		     [&pop]( const typename sugarpop_t::mutation_t & __m ) { pop.mut_lookup.insert(__m.pos); } );

    }

   template<typename sugarpop_t,
	     typename reader_t,
	     typename diploid_reader_t = diploidIOplaceholder>
    inline typename std::enable_if<std::is_same<typename sugarpop_t::popmodel_t,sugar::MULTILOCPOP_TAG>::value,result_type>::type
    operator()( sugarpop_t & pop,
		const serialize & s,
		const reader_t & rt,
		const diploid_reader_t & dr = diploid_reader_t()) const
    {
      pop.clear();
      //Step 0: read N
      s.buffer.read( reinterpret_cast<char*>(&pop.N),sizeof(unsigned) );
      KTfwd::read_binary_pop_mloc( &pop.gametes,&pop.mutations,&pop.diploids,rt,s.buffer,dr );
      unsigned temp;
      s.buffer.read( reinterpret_cast<char*>(&temp),sizeof(unsigned) );
      for( unsigned m=0;m<temp ;++m )
    	{
    	  typename reader_t::result_type mm = rt(s.buffer);
    	  pop.fixations.emplace_back( std::move(mm) );
    	}
      pop.fixation_times.resize(temp);
      if(temp)
	{
	  s.buffer.read( reinterpret_cast<char*>(&pop.fixation_times[0]), temp*sizeof(unsigned) );
	}
      s.buffer.seekg(0);
      //Finally, fill the lookup table:
      std::for_each( pop.mutations.begin(), pop.mutations.end(),
    		     [&pop]( const typename sugarpop_t::mutation_t & __m ) { pop.mut_lookup.insert(__m.pos); } );

    }

    /*!
      \brief Overload for metapopulation simulations
     */
    template<typename sugarpop_t,
	     typename reader_t,
	     typename diploid_reader_t = diploidIOplaceholder>
    inline typename std::enable_if<std::is_same<typename sugarpop_t::popmodel_t,sugar::METAPOP_TAG>::value,result_type>::type
    operator()( sugarpop_t & pop,
		const serialize & s,
		const reader_t & rt,
		const diploid_reader_t & dr = diploid_reader_t()) const
    {
      pop.clear();
      //Step 0: read N
      unsigned numNs;
      s.buffer.read( reinterpret_cast<char*>(&numNs),sizeof(unsigned) );
      pop.Ns.resize(numNs);
      s.buffer.read( reinterpret_cast<char*>(&pop.Ns[0]),numNs*sizeof(unsigned) );
      //Step 1: write the mutations, diploids, gametes to the stream
      KTfwd::read_binary_metapop( &pop.gametes,&pop.mutations,&pop.diploids,rt,s.buffer,dr );
      unsigned temp;
      s.buffer.read( reinterpret_cast<char*>(&temp),sizeof(unsigned) );
      for( unsigned m=0;m<temp ;++m )
    	{
    	  typename reader_t::result_type mm = rt(s.buffer);
    	  pop.fixations.emplace_back( std::move(mm) );
    	}
      pop.fixation_times.resize(temp);
      if(temp)
	{
	  s.buffer.read( reinterpret_cast<char*>(&pop.fixation_times[0]), temp*sizeof(unsigned) );
	}
      s.buffer.seekg(0);
      //Finally, fill the lookup table:
      std::for_each( pop.mutations.begin(), pop.mutations.end(),
    		     [&pop]( const typename sugarpop_t::mutation_t & __m ) { pop.mut_lookup.insert(__m.pos); } );
    }
  };

  /*!
    Write a population to a gzFile in a binary format
  */
  struct gzserialize
  {
    using result_type = long long;
    /*
      \brief Call operator
      \note gzout must already be opened, and with a mode involving 'b'
     */
    template<typename sugarpop_t,
	     typename writer_t,
	     typename diploid_writer_t = diploidIOplaceholder>
    inline long long operator()( gzFile gzout,
				 const sugarpop_t & pop,
				 const writer_t & wt,
				 const diploid_writer_t & dw = diploid_writer_t() ) const
    {
      serialize s;
      s(pop,wt,dw);
      return gzwrite(gzout,s.buffer.str().c_str(),s.buffer.str().size());
    }
  };

  /*!
    Read a population from a gzFile in binary format
   */
  struct gzdeserialize
  {
    using result_type = void;
    /*!
      \brief Call operator
      \note gzin must be opened for reading in binary mode
    */
    template<typename sugarpop_t,
	     typename reader_t,
	     typename diploid_reader_t = diploidIOplaceholder>
    inline typename std::enable_if<std::is_same<typename sugarpop_t::popmodel_t,sugar::SINGLEPOP_TAG>::value,
				   result_type>::type operator()( gzFile gzin,
								  sugarpop_t & pop,
								  const reader_t & rt,
								  const diploid_reader_t & dr = diploid_reader_t()) const
    {
      pop.clear();
      //Step 0: read N
      gzread( gzin, reinterpret_cast<char*>(&pop.N),sizeof(unsigned) );
      KTfwd::read_binary_pop( &pop.gametes,&pop.mutations,&pop.diploids,rt,gzin,dr );
      unsigned temp;
      gzread( gzin, reinterpret_cast<char*>(&temp),sizeof(unsigned) );
      for( unsigned m=0;m<temp ;++m )
    	{
    	  typename reader_t::result_type mm = rt(gzin);
    	  pop.fixations.emplace_back( std::move(mm) );
    	}
      pop.fixation_times.resize(temp);
      if(temp)
	{
	  gzread( gzin, reinterpret_cast<char*>(&pop.fixation_times[0]), temp*sizeof(unsigned) );
	}
      //Finally, fill the lookup table:
      std::for_each( pop.mutations.begin(), pop.mutations.end(),
    		     [&pop]( const typename sugarpop_t::mutation_t & __m ) { pop.mut_lookup.insert(__m.pos); } );
    }

    template<typename sugarpop_t,
	     typename reader_t,
	     typename diploid_reader_t = diploidIOplaceholder>
    inline typename std::enable_if<std::is_same<typename sugarpop_t::popmodel_t,sugar::MULTILOCPOP_TAG>::value,
				   result_type>::type operator()( gzFile gzin,
								  sugarpop_t & pop,
								  const reader_t & rt,
								  const diploid_reader_t & dr = diploid_reader_t()) const
    {
      pop.clear();
      //Step 0: read N
      gzread( gzin, reinterpret_cast<char*>(&pop.N),sizeof(unsigned) );
      KTfwd::read_binary_pop_mloc( &pop.gametes,&pop.mutations,&pop.diploids,rt,gzin,dr );
      unsigned temp;
      gzread( gzin, reinterpret_cast<char*>(&temp),sizeof(unsigned) );
      for( unsigned m=0;m<temp ;++m )
    	{
    	  typename reader_t::result_type mm = rt(gzin);
    	  pop.fixations.emplace_back( std::move(mm) );
    	}
      pop.fixation_times.resize(temp);
      if(temp)
	{
	  gzread( gzin, reinterpret_cast<char*>(&pop.fixation_times[0]), temp*sizeof(unsigned) );
	}
      //Finally, fill the lookup table:
      std::for_each( pop.mutations.begin(), pop.mutations.end(),
    		     [&pop]( const typename sugarpop_t::mutation_t & __m ) { pop.mut_lookup.insert(__m.pos); } );
    }

    /*!
      \brief Call operator
      \note gzin must be opened for reading in binary mode
    */
    template<typename sugarpop_t,
	     typename reader_t,
	     typename diploid_reader_t = diploidIOplaceholder>
    inline typename std::enable_if<std::is_same<typename sugarpop_t::popmodel_t,sugar::METAPOP_TAG>::value,result_type>::type
    operator()( gzFile gzin,
		sugarpop_t & pop,
		const reader_t & rt,
		const diploid_reader_t & dr = diploid_reader_t()) const
    {
      pop.clear();
      //Step 0: read N
      unsigned numNs;
      gzread( gzin, reinterpret_cast<char*>(&numNs),sizeof(unsigned) );
      pop.Ns.resize(numNs);
      gzread( gzin, reinterpret_cast<char*>(&pop.Ns[0]),numNs*sizeof(unsigned) );
      //Step 1: write the mutations, diploids, gametes to the stream
      KTfwd::read_binary_metapop( &pop.gametes,&pop.mutations,&pop.diploids,rt,gzin,dr );
      unsigned temp;
      gzread( gzin, reinterpret_cast<char*>(&temp),sizeof(unsigned) );
      for( unsigned m=0;m<temp ;++m )
    	{
    	  typename reader_t::result_type mm = rt(gzin);
    	  pop.fixations.emplace_back( std::move(mm) );
    	}
      pop.fixation_times.resize(temp);
      if(temp)
	{
	  gzread( gzin, reinterpret_cast<char*>(&pop.fixation_times[0]), temp*sizeof(unsigned) );
	}
      //Finally, fill the lookup table:
      std::for_each( pop.mutations.begin(), pop.mutations.end(),
    		     [&pop]( const typename sugarpop_t::mutation_t & __m ) { pop.mut_lookup.insert(__m.pos); } );
    }
  };
}

#endif
