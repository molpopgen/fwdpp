/*
Essentially a re-implementation of diploid_ind.cc that is capable of being run in python.

In addition to the usual fwdpp depdencies, we need boost.python.

To compile:
g++ -fPIC -Wall -W -O3 `python-config --includes` -std=c++11 -c fwdpp_boost_python.cc
g++ -std=c++11 -shared -o fwdpp_boost_python.so fwdpp_boost_python.o -lboost_python -lboost_system  -lpython -lgsl -lgslcblas

To run:
python test_boost_python.py
*/

//include main fwdpp library
#include <fwdpp/diploid.hh>
//Tell the fwdpp "sugar" layer to use boost containers (no harm here, as it is likely installed b/c we need boost python..)
#define FWDPP_SUGAR_USE_BOOST
//Include the necessary "sugar" components
#include <fwdpp/sugar/singlepop.hpp>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <boost/python.hpp>

using namespace boost::python;

/*
  We will use a gsl_rng_mt19937 as our RNG.
  This type is implicitly convertible to gsl_rng *,
  and auto-handles the gsl_rng_free steps, etc.
 */
using GSLrng = KTfwd::GSLrng_t<KTfwd::GSL_RNG_MT19937>;

/*
  boost.python requires that any C++ classes exposed to python
  are copy-constructible.

  As far as fwdpp is concerned, this requirement means that 
  we must use the serializable version of the fwdpp sugar libary.

  In order to use those object, we must define objects that will write
  and read mutation types.
 */
struct mwriter
{
  using result_type = void;
  inline result_type operator()( const KTfwd::mutation & m,
				 std::ostream & o ) const
  {
    o.write( reinterpret_cast<const char *>(&m.pos),sizeof(double));
    o.write( reinterpret_cast<const char *>(&m.n),sizeof(unsigned));
    o.write( reinterpret_cast<const char *>(&m.neutral),sizeof(bool));
    o.write( reinterpret_cast<const char *>(&m.s),sizeof(double));
    o.write( reinterpret_cast<const char *>(&m.h),sizeof(double));
  }
};

struct mreader
{
  using result_type = KTfwd::mutation;
  inline result_type operator()( std::istream & in ) const
  {
    double pos,s,h;
    unsigned n;
    bool neutral;
    in.read( reinterpret_cast< char *>(&pos),sizeof(double));
    in.read( reinterpret_cast< char *>(&n),sizeof(unsigned));
    in.read( reinterpret_cast< char *>(&neutral),sizeof(bool));
    in.read( reinterpret_cast< char *>(&s),sizeof(double));
    in.read( reinterpret_cast< char *>(&h),sizeof(double));
    return result_type(pos,s,n,h);
  }
};

/*
  With the above definitions, we can not define what our population looks like.
 */

using poptype = KTfwd::singlepop_serialized<KTfwd::mutation,mwriter,mreader>;

//mutation function is infinitely-many sites, all neutral variants
poptype::mtype neutral_mutations_inf_sites(gsl_rng * r,
					   poptype::mlist_t * mutations,
					   poptype::lookup_table_t * lookup)
{
  double pos = gsl_rng_uniform(r);
  while( lookup->find(pos) != lookup->end() ) //make sure it doesn't exist in the population
    { 
      pos = gsl_rng_uniform(r);  //if it does, generate a new one
    }
  lookup->insert(pos);
  return poptype::mtype(pos,0.,1,true);
}

//Evolve the population for some amount of time with mutation and recombination
poptype evolve( GSLrng & rng,
		const unsigned & N,
		const unsigned & generations,
		const double & mu,
		const double & recrate)
{
  poptype pop(N);
  std::function<double(void)> recmap = std::bind(gsl_rng_uniform,rng); //uniform crossover map
  for( unsigned generation = 0 ; generation < generations ; ++generation )
    {
      double wbar = KTfwd::sample_diploid(rng,
					  &pop.gametes,
					  &pop.diploids,
					  &pop.mutations,
					  N,
					  mu,
					  std::bind(neutral_mutations_inf_sites,rng,std::placeholders::_1,&pop.mut_lookup),
					  std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,
						    &pop.gametes,
						    recrate, 
						    rng,
						    recmap),
					  std::bind(KTfwd::insert_at_end<poptype::mtype,poptype::mlist_t>,std::placeholders::_1,std::placeholders::_2),
					  std::bind(KTfwd::insert_at_end<poptype::gtype,poptype::glist_t>,std::placeholders::_1,std::placeholders::_2),
					  std::bind(KTfwd::multiplicative_diploid(),std::placeholders::_1,std::placeholders::_2,2.),
					  std::bind(KTfwd::mutation_remover(),std::placeholders::_1,0,2*pop.N));
      KTfwd::remove_fixed_lost(&pop.mutations,&pop.fixations,&pop.fixation_times,&pop.mut_lookup,generation,2*pop.N);
    }
  return pop;
}

//Calculate the site-frequency spectrum for a sample
boost::python::list sfs(GSLrng & rng,const poptype & pop,const unsigned & nsam)
{
  std::map<double,unsigned> mutfreqs;
  unsigned twoN = 2*pop.N;

  for( unsigned i = 0 ; i < nsam ; ++i )
    {
      //pick a random chrom (w/replacement...)
      unsigned chrom = unsigned(gsl_ran_flat(rng,0.,double(twoN)));
      //get pointer to that chrom from the individual
      auto gamete = (chrom%2==0.) ? pop.diploids[chrom/2].first : pop.diploids[chrom/2].second;
      //In this example, there are only neutral mutations, so that's what we'll iterate over
      for( auto m = gamete->mutations.begin() ; m != gamete->mutations.end() ; ++m )
	{
	  auto pos_itr = mutfreqs.find( (*m)->pos );
	  if( pos_itr == mutfreqs.end() )
	    {
	      mutfreqs.insert(std::make_pair((*m)->pos,1));
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

