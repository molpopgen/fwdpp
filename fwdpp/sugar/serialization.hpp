#ifndef __FWDPP_SUGAR_SINGLEPOP_SERIALIZATON_HPP__
#define __FWDPP_SUGAR_SINGLEPOP_SERIALIZATON_HPP__

#include <sstream>
#include <algorithm>
#include <type_traits>
#include <fwdpp/IO.hpp>
#include <fwdpp/sugar/poptypes/tags.hpp>

namespace KTfwd
{
  struct serialize
  {
    using result_type = void;
    mutable std::ostringstream buffer;
    template<typename sugarpop_t,
	     typename writer_t>
    inline typename std::enable_if<std::is_same<typename sugarpop_t::popmodel_t,sugar::SINGLEPOP_TAG>::value,result_type>::type
    operator()( const sugarpop_t & pop,
		const writer_t & wt ) const
    {
      buffer.str(std::string());
      buffer.write( reinterpret_cast<const char*>(&pop.N), sizeof(unsigned) );
      write_binary_pop(&pop.gametes,&pop.mutations,&pop.diploids,wt,buffer);
      //Step 2: output fixations 
      unsigned temp = pop.fixations.size();
      buffer.write( reinterpret_cast<char*>(&temp), sizeof(unsigned) );
      std::for_each( pop.fixations.begin(), pop.fixations.end(),
		     std::bind(wt,std::placeholders::_1,std::ref(buffer)) );
      //Step 3:the fixation times
      buffer.write( reinterpret_cast<const char *>(&pop.fixation_times[0]), pop.fixation_times.size()*sizeof(unsigned) );
    }

    template<typename sugarpop_t,
	     typename writer_t>
    inline typename std::enable_if<std::is_same<typename sugarpop_t::popmodel_t,sugar::METAPOP_TAG>::value,result_type>::type
    operator()( const sugarpop_t & pop,
		const writer_t & wt ) const
    {
      buffer.str(std::string());
      unsigned npops = pop.Ns.size();
      buffer.write(reinterpret_cast<char *>(&npops),sizeof(unsigned));
      buffer.write( reinterpret_cast<const char*>(&pop.Ns[0]), npops*sizeof(unsigned) );
      write_binary_pop(&pop.gametes,&pop.mutations,&pop.diploids,wt,buffer);
      //Step 2: output fixations 
      unsigned temp = pop.fixations.size();
      buffer.write( reinterpret_cast<char*>(&temp), sizeof(unsigned) );
      std::for_each( pop.fixations.begin(), pop.fixations.end(),
		     std::bind(wt,std::placeholders::_1,std::ref(buffer)) );
      //Step 3:the fixation times
      buffer.write( reinterpret_cast<const char *>(&pop.fixation_times[0]), pop.fixation_times.size()*sizeof(unsigned) );
    }
			  
  };

  struct deserialize
  {
    using result_type = void;
    template<typename sugarpop_t,
	     typename reader_t>
    inline typename std::enable_if<std::is_same<typename sugarpop_t::popmodel_t,sugar::SINGLEPOP_TAG>::value,result_type>::type
    operator()( sugarpop_t & pop,
		const serialize & s,
		const reader_t & rt ) const
    {
      pop.clear();
      std::istringstream i(s.buffer.str());
      //Step 0: read N
      i.read( reinterpret_cast<char*>(&pop.N),sizeof(unsigned) );
      //Step 1: write the mutations, diploids, gametes to the stream
      KTfwd::read_binary_pop( &pop.gametes,&pop.mutations,&pop.diploids,rt,i );
      unsigned temp;
      i.read( reinterpret_cast<char*>(&temp),sizeof(unsigned) );
      for( unsigned m=0;m<temp ;++m )
    	{
    	  typename reader_t::result_type mm = rt(i);
    	  pop.fixations.emplace_back( std::move(mm) );
    	}
      pop.fixation_times.resize(temp);
      i.read( reinterpret_cast<char*>(&pop.fixation_times[0]), temp*sizeof(unsigned) );

      //Finally, fill the lookup table:
      std::for_each( pop.mutations.begin(), pop.mutations.end(),
    		     [&pop]( const typename sugarpop_t::mutation_t & __m ) { pop.mut_lookup.insert(__m.pos); } );
    }
  };
}

#endif
