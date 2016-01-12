/*
Essentially a re-implementation of diploid_ind.cc that is capable of being run in python.

In addition to the usual fwdpp depdencies, we need boost.python.

To compile:
g++ -fPIC -Wall -W -O3 -I. `python-config --includes` -std=c++11 -c fwdpp_boost_python.cc
g++ -std=c++11 -shared -o fwdpp_boost_python.so fwdpp_boost_python.o -lboost_python -lboost_system  -lpython -lgsl -lgslcblas

To run:
python test_boost_python.py
*/
#include <python_common.hpp>
#include <boost/python.hpp>

using namespace boost::python;

//Calculate the site-frequency spectrum for a sample
boost::python::list sfs(GSLrng & rng,const poptype & pop,const unsigned & nsam)
{
  std::map<double,unsigned> mutfreqs;
  unsigned twoN = 2*pop.N;

  for( unsigned i = 0 ; i < nsam ; ++i )
    {
      //pick a random chrom (w/replacement...)
      unsigned chrom = unsigned(gsl_ran_flat(rng.get(),0.,double(twoN)));
      //get reference to that chrom from the individual
      auto & gamete = (chrom%2==0.) ? pop.gametes[pop.diploids[chrom/2].first] : pop.gametes[pop.diploids[chrom/2].second];
      //In this example, there are only neutral mutations, so that's what we'll iterate over
      for(auto & m : gamete.mutations)
	{
	  auto pos = pop.mutations[m].pos;
	  auto pos_itr = mutfreqs.find( pos );
	  if( pos_itr == mutfreqs.end() )
	    {
	      mutfreqs.insert(std::make_pair(pos,1));
	    }
	  else
	    {
	      pos_itr->second++;
	    }
	}
    }
  //Now, fill in the SFS, omitting positions that are fixed in the sample
  std::vector<unsigned> __rv(nsam-1,0u);
  for( const auto & __x : mutfreqs )
    {
      if (__x.second < nsam) __rv[__x.second-1]++;
    }
  boost::python::list rv;
  for( const auto & __x : __rv ) rv.append(__x);
  return rv;
}

//Now, we can expose the stuff to python
BOOST_PYTHON_MODULE(fwdpp_boost_python)
{
  //Expose the type based on fwdpp's "sugar" layer
  class_<poptype>("poptype",init<unsigned>())
    .def("clear",&poptype::clear)
    ;
  //Expose the GSL wrapper
  class_<GSLrng>("GSLrng",init<unsigned>())
    ;
  //Expose the function to run the model
  def("evolve",evolve);
  //And one to get the sfs of a sample
  def("sfs",sfs);
}

