/*! \include migsel_ind.cc
  Two constant-size populations with selection and inbreeding.

  The selection coefficient, s, is treated as -s in deme 2, just for fun.

  Writes the metapop + an "ms"-type sample in binary format to an output file.
 */
#include <config.h>
#include <fwdpp/diploid.hh>

#include <Sequence/SimData.hpp>
#include <Sequence/SimDataIO.hpp> //for writing & reading SimData objects in binary format
#include <Sequence/FST.hpp>
#include <numeric>
#include <cmath>
#include <functional>
#include <cassert>
#include <iomanip>
#include <fstream>
//the type of mutation
typedef KTfwd::mutation mtype;
#include <common_ind.hpp>

using namespace std;
using namespace KTfwd;
using namespace Sequence;

//function object to write mutation data in binary format
struct mwriter
{
  typedef void result_type;
  result_type operator()( const mtype & m, std::ostream & buffer ) const
  {
    unsigned u = m.n;
    buffer.write( reinterpret_cast< char * >(&u),sizeof(unsigned) );
    bool b = m.neutral;
    buffer.write( reinterpret_cast< char * >(&b),sizeof(bool) );
    double d = m.pos;
    buffer.write( reinterpret_cast< char * >(&d),sizeof(double) );
    d = m.s;
    buffer.write( reinterpret_cast< char * >(&d),sizeof(double) );
    d = m.h;
    buffer.write( reinterpret_cast< char * >(&d),sizeof(double) );
  }
};

//function object to read mutation data in binary format
struct mreader
{
  typedef mtype result_type;
  result_type operator()( std::istream & in ) const
  {
    unsigned n;
    in.read( reinterpret_cast< char * >(&n),sizeof(unsigned) );
    bool neut;
    in.read( reinterpret_cast< char * >(&neut),sizeof(bool) );
    double pos;
    in.read( reinterpret_cast< char * >(&pos),sizeof(double) );
    double s;
    in.read( reinterpret_cast< char * >(&s),sizeof(double) );
    double h;
    in.read( reinterpret_cast< char * >(&h),sizeof(double) );
    return result_type(pos,n,s,h);
  }
};

mutation neutral_mutations_selection(gsl_rng * r,mlist * mutations,
					    const double & mu_neutral, const double & mu_selected,
					    const double & s, const double & h,
					    lookup_table_type * lookup)
{
  /*
    Choose a mutation position [0,1) that does not currently exist in the population.
  */
  double pos = gsl_rng_uniform(r);
  while( lookup->find(pos) != lookup->end() )
    {
      pos = gsl_rng_uniform(r);
    }
  lookup->insert(pos);
  assert(std::find_if(mutations->begin(),mutations->end(),std::bind(mutation_at_pos(),std::placeholders::_1,pos)) == mutations->end());
  //Choose mutation class (neutral or selected) proportional to underlying mutation rates
  bool neutral_mut = ( gsl_rng_uniform(r) <= mu_neutral/(mu_neutral+mu_selected) ) ? true : false;
  //Return a mutation of the correct type.  Neutral means s = h = 0, selected means s=s and h=h.
  return (neutral_mut == true) ? mutation(pos,0.,1,0.) : mutation(pos,s,1,h);
}

size_t migpop(const size_t & source_pop, gsl_rng * r, const double & mig_prob)
{
  if( gsl_rng_uniform(r) <= mig_prob )
    {
      return ! source_pop;
    }
  return source_pop;
}

SimData merge( const std::vector<std::pair<double,std::string> > & sample1,
			 const std::vector<std::pair<double,std::string> > & sample2 ,
			 const unsigned & nsam);

//fitness model details -- s will be treated as -s in population 2
struct multiplicative_fitness_updater_hom_minus
{
  typedef void result_type;
  template<typename iterator_type>
  inline void operator()(double & fitness, const iterator_type & m1,const double & scaling = 2.) const
  {
    fitness *= ( 1. - scaling*m1->s );
  }
};

struct multiplicative_fitness_updater_het_minus
{
  typedef void result_type;
  template<typename iterator_type>
  inline void operator()(double & fitness, const iterator_type & m1) const
  {
    fitness *= ( 1. - m1->s*m1->h );
  }
};

struct multiplicative_diploid_minus
{
  typedef double result_type;
  template< typename iterator_type>
  inline double operator()(const iterator_type & g1, const iterator_type & g2, 
			   const double scaling = 1.) const
  {
    return site_dependent_fitness()(g1,g2,
				    std::bind(multiplicative_fitness_updater_hom_minus(),std::placeholders::_1,std::placeholders::_2,scaling),
				    std::bind(multiplicative_fitness_updater_het_minus(),std::placeholders::_1,std::placeholders::_2),
				    1.);
  }
};

int main( int argc, char ** argv )
{
  if (argc != 14)
    {
      std::cerr << "Too few arguments.\n"
		<< "Usage: migsel_ind N theta_neutral theta_deleterious rho M s h f1 f2 ngens n outfilename seed\n";
      exit(10);
    }
  int argn=1;
  const unsigned N = atoi(argv[argn++]);
  const double theta_neut = atof(argv[argn++]);
  const double theta_del = atof(argv[argn++]);
  const double rho = atof(argv[argn++]);
  const double M = atof(argv[argn++]);
  const double s = atof(argv[argn++]);
  const double h = atof(argv[argn++]);
  const double f1 = atof(argv[argn++]);
  const double f2 = atof(argv[argn++]);
  const unsigned ngens = atoi(argv[argn++]);
  const unsigned n = atoi(argv[argn++]);
  const char * outfilename = argv[argn++];
  const unsigned seed = atoi(argv[argn++]);

  const double mu_neutral = theta_neut/double(4*N);
  const double mu_del = theta_del/double(4*N);
  const double littler = rho/double(4*N);
  const double m = M/double(4*N);

  gsl_rng * r =  gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(r,seed);

  mlist mutations; 
  std::vector<mtype> fixations;
  std::vector<unsigned> fixation_times;

  lookup_table_type lookup;
#if defined(HAVE_BOOST_VECTOR) && defined(HAVE_BOOST_LIST) && defined(HAVE_BOOST_UNORDERED_SET) && defined(HAVE_BOOST_POOL_ALLOC) && defined(HAVE_BOOST_HASH) && !defined(USE_STANDARD_CONTAINERS)
  boost::container::vector< glist > metapop(2, glist(1,gtype(2*N)));
  typedef boost::container::vector< std::pair< glist::iterator, glist::iterator > > diploid_bucket;
  boost::container::vector< diploid_bucket > diploids;
  boost::container::vector<unsigned> Ns(2,N);
  boost::container::vector<double> fs;
#else
  std::vector< glist > metapop(2, glist(1,gtype(2*N)));
  typedef std::vector< std::pair< glist::iterator, glist::iterator > > diploid_bucket;
  std::vector< diploid_bucket > diploids;
  std::vector<unsigned> Ns(2,N);
  std::vector<double> fs;
#endif


  for ( auto i = metapop.begin() ;
	i != metapop.end() ; ++i )
    {
      diploids.push_back ( diploid_bucket(N,std::make_pair(i->begin(),i->begin())) );
    }

  fs.push_back(f1);
  fs.push_back(f2);

  //create a vector of fitness functions for each population
  std::vector<std::function<double (glist::const_iterator,
				      glist::const_iterator)> > vbf;
  vbf.push_back(std::bind(multiplicative_diploid(),std::placeholders::_1,std::placeholders::_2,2.));
  vbf.push_back(std::bind(multiplicative_diploid_minus(),std::placeholders::_1,std::placeholders::_2,2.));

  //recombination map is uniform[0,1)
  std::function<double(void)> recmap = std::bind(gsl_rng_uniform,r);

  for( unsigned generation = 0 ; generation < ngens ; ++generation )
    {
      std::vector<double> wbars = sample_diploid(r,
						 &metapop,
						 &diploids,
						 &mutations,
						 &Ns[0],
						 mu_neutral + mu_del,
						 std::bind(neutral_mutations_selection,r,std::placeholders::_1,mu_neutral,mu_del,s,h,&lookup),
						 std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,
							     littler,
							     r,
							     recmap),
						 std::bind(insert_at_end<mtype,mlist>,std::placeholders::_1,std::placeholders::_2),
						 std::bind(insert_at_end<gtype,glist>,std::placeholders::_1,std::placeholders::_2),
						 vbf,
						 //4*N b/c it needs to be fixed in the metapopulation
						 std::bind(mutation_remover(),std::placeholders::_1,0,4*N),
						 std::bind(migpop,std::placeholders::_1,r,m),
						 &fs[0]);
      //4*N b/c it needs to be fixed in the metapopulation
      remove_fixed_lost(&mutations,&fixations,&fixation_times,&lookup,generation,4*N);
    }
  
  std::pair< std::vector<std::pair<double,std::string> >,
	     std::vector<std::pair<double,std::string> > > spop1 = ms_sample_separate(r,&diploids[0],n);

  std::pair< std::vector<std::pair<double,std::string> >,
	     std::vector<std::pair<double,std::string> > > spop2 = ms_sample_separate(r,&diploids[1],n);

  SimData neutral = merge( spop1.first,spop2.first,n );
  SimData selected = merge( spop1.second,spop2.second,n );

  std::ofstream outstream(outfilename);

  //Write the metapop in binary format to outstream
  KTfwd::write_binary_metapop(&metapop,&mutations,&diploids,
   			      std::bind(mwriter(),std::placeholders::_1,std::placeholders::_2),
			      outstream);

  //Write the "ms" blocks
  Sequence::write_SimData_binary(outstream,neutral);
  Sequence::write_SimData_binary(outstream,selected);
  
  outstream.close();

  //Now, read it all back in, for fun.
#if defined(HAVE_BOOST_VECTOR) && defined(HAVE_BOOST_LIST) && defined(HAVE_BOOST_UNORDERED_SET) && defined(HAVE_BOOST_POOL_ALLOC) && defined(HAVE_BOOST_HASH) && !defined(USE_STANDARD_CONTAINERS)
  boost::container::vector< glist > metapop2;
  boost::container::vector< diploid_bucket > diploids2;
#else
  std::vector< glist > metapop2;
  std::vector< diploid_bucket > diploids2;
#endif
  mlist mutations2;
  Sequence::SimData neutral2,selected2;

  ifstream in(outfilename);
  
  KTfwd::read_binary_metapop(&metapop2,&mutations2,&diploids2,
			     std::bind(mreader(),std::placeholders::_1),
			     in);
  
  assert( metapop2.size() == metapop.size() );
  assert( mutations2.size() == mutations.size() );
  assert( diploids2.size() == diploids.size() );
  
  neutral2 = Sequence::read_SimData_binary(in);
  selected2 = Sequence::read_SimData_binary(in);
  
  in.close();
  
  std::cerr << (neutral == neutral2) << ' ' << (selected == selected2) << '\n';
  /*
    At this point, you could go through each deme and each diploid and make 
    sure that all is cool.  However, if we weren't reading and
    writing the metapop correctly, there's no way we'd be able
    to write and then read in the ms blocks correctly, as we'd have
    run into some binary gibberish along the way.
  */

  //For fun, we'll calculate some basic pop subdivision stats
  unsigned config[2] = {n,n};
  if(!neutral.empty())
    {
      Sequence::FST fst_neut(&neutral,2,config);
      std::pair< std::set<double>,std::set<double> > pneut = fst_neut.Private(0,1);
      std::cout << fst_neut.HSM() << '\t'
		<< fst_neut.shared(0,1).size() << '\t'
		<< pneut.first.size() << '\t'
		<< pneut.second.size() << '\t';
    }
  else
    {
      std::cout << "NA\t0\t0\t0\t0\t";
    }
  if(!selected.empty())
    {
      Sequence::FST fst_sel(&selected,2,config);
      std::pair< std::set<double>,std::set<double> > psel = fst_sel.Private(0,1);
      std::cout << fst_sel.HSM() << '\t'
		<< fst_sel.shared(0,1).size() << '\t'
		<< psel.first.size() << '\t'
		<< psel.second.size() << '\n';    
    }
  else
    {
      std::cout << "NA\t0\t0\t0\t0\n";
    }
}

SimData merge( const std::vector<std::pair<double,std::string> > & sample1,
			 const std::vector<std::pair<double,std::string> > & sample2 ,
			 const unsigned & nsam)
{
  std::map<double, std::string> temp;

  for( unsigned i=0;i<sample1.size();++i)
    {
      temp[sample1[i].first] = std::string(sample1[i].second + std::string(nsam,'0'));
    }

  for( unsigned i=0;i<sample2.size();++i)
    {
      std::map<double,std::string>::iterator itr = temp.find(sample2[i].first);
      if( itr == temp.end() )
	{
	  temp[sample2[i].first] = std::string( std::string(nsam,'0') + sample2[i].second );
	}
      else
	{
	  std::copy( sample2[i].second.begin(),sample2[i].second.end(),itr->second.begin()+nsam );
	}
    }
  std::vector<std::pair<double,std::string> > rv( temp.begin(),temp.end() );
  std::sort(rv.begin(),rv.end(),
		[](std::pair<double,std::string> lhs,
		   std::pair<double,std::string> rhs) { return lhs.first < rhs.first; });
  return SimData(rv.begin(),rv.end());
}
