/*! \include migsel_split_ind_list.cc
  Evolve a constant-N pop, split it off, then evolve both with migration.

  The selection coefficient, s, is treated as -s in deme 2, just for fun.

  Output is a "ms" block for neutral and non-neutral variants.  Ancestral
  population first.
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
#include <iostream>
#include <fwdpp/sugar/infsites.hpp>
#include <fwdpp/sugar/metapop/metapop.hpp>

#if defined(HAVE_BOOST_VECTOR) && defined(HAVE_BOOST_LIST) && defined(HAVE_BOOST_UNORDERED_SET) && defined(HAVE_BOOST_POOL_ALLOC) && defined(HAVE_BOOST_HASH) && !defined(USE_STANDARD_CONTAINERS)
#include <boost/container/list.hpp>
#include <boost/container/vector.hpp>
#include <boost/pool/pool_alloc.hpp>
#include <boost/unordered_set.hpp>
#include <boost/functional/hash.hpp>

template<typename mtype> using metapop_mlist_t = boost::container::list<mtype,boost::fast_pool_allocator<mtype> >;
template<typename mtype> using metapop_gamete_t = KTfwd::gamete_base<mtype,metapop_mlist_t<mtype>>;
template<typename mtype> using metapop_glist_t = boost::container::list<metapop_gamete_t<mtype>, boost::fast_pool_allocator<metapop_gamete_t<mtype>>>;
template<typename mtype> using metapop_dipvector_t = boost::container::vector<std::pair<typename metapop_glist_t<mtype>::iterator,
											typename metapop_glist_t<mtype>::iterator> >;

template<typename mtype> using metapop_t = KTfwd::sugar::metapop<mtype,
							       metapop_mlist_t<mtype>,
							       metapop_glist_t<mtype>,
							       metapop_dipvector_t<mtype>,
							       boost::container::list<metapop_glist_t<mtype>>,
							       boost::container::vector<metapop_dipvector_t<mtype>>,
							       boost::container::vector<mtype>,
							       boost::container::vector<unsigned>,
							       boost::unordered_set<double,boost::hash<double>,KTfwd::equal_eps>>;

#else
#include <list>
#include <vector>
#include <unordered_set>

template<typename mtype> using metapop_mlist_t = std::list<mtype>;
template<typename mtype> using metapop_gamete_t = KTfwd::gamete_base<mtype,metapop_mlist_t<mtype>>;
template<typename mtype> using metapop_glist_t = std::list<metapop_gamete_t<mtype> >;
template<typename mtype> using metapop_dipvector_t = std::vector<std::pair<typename metapop_glist_t<mtype>::iterator,
									   typename metapop_glist_t<mtype>::iterator> >;

template<typename mtype> using metapop_t = KTfwd::sugar::metapop<mtype,
							       metapop_mlist_t<mtype>,
							       metapop_glist_t<mtype>,
							       metapop_dipvector_t<mtype>,
							       std::list<metapop_glist_t<mtype>>,
							       std::vector<metapop_dipvector_t<mtype>>,
							       std::vector<mtype>,
							       std::vector<unsigned>,
							       std::unordered_set<double,std::hash<double>,KTfwd::equal_eps>>;
#endif


using mtype = KTfwd::mutation;

using namespace std;
using namespace KTfwd;
using namespace Sequence;

using poptype = metapop_t<mtype>;
using glist = poptype::glist_t;
using mlist = poptype::mlist_t;
using gtype = poptype::gamete_t;
using diploid_bucket = poptype::dipvector_t;
using diploid_bucket_vec = poptype::vdipvector_t;
using glist_list = poptype::vglist_t;

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
    return std::max(0.,site_dependent_fitness()(g1,g2,
						std::bind(multiplicative_fitness_updater_hom_minus(),std::placeholders::_1,std::placeholders::_2,scaling),
						std::bind(multiplicative_fitness_updater_het_minus(),std::placeholders::_1,std::placeholders::_2),
						1.));
  }
};

void duplicate_pop(glist_list & metapop, diploid_bucket_vec & diploids, mlist & mutations)
{
  diploids.push_back( diploid_bucket(diploids[0]) );
  metapop.insert( metapop.end(), glist(*metapop.begin()) );
  auto pptr1 = metapop.begin(),pptr2=metapop.begin();
  ++pptr2;
  for (auto dptr=diploids[0].begin(), newdptr=diploids[1].begin(); dptr!=diploids[0].end(); ++dptr, ++newdptr)
    {
      auto pos1=std::distance(pptr1->begin(), dptr->first);
      auto pos2=std::distance(pptr1->begin(), dptr->second);
      auto newgptr1=pptr2->begin();
      std::advance(newgptr1, pos1);
      auto newgptr2=pptr2->begin();
      std::advance(newgptr2, pos2);
      newdptr->first=newgptr1;
      newdptr->second=newgptr2;
    }
  for (auto mptr=mutations.begin(); mptr!=mutations.end(); ++mptr)
    {
      mptr->n+=mptr->n;
    }
}

int main( int argc, char ** argv )
{
  if (argc != 14)
    {
      std::cerr << "Too few arguments.\n"
		<< "Usage: " << argv[0] << " N theta_neutral theta_deleterious rho M s h f1 f2 ngens ngens2 nsam seed\n";
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
  const unsigned ngens2 = atoi(argv[argn++]);
  const unsigned n = atoi(argv[argn++]);
  const unsigned seed = atoi(argv[argn++]);

  if(!ngens2) 
    {
      std::cerr << "Error: ngens2 must be > 0\n";
      std::exit(EXIT_FAILURE);
    }

  const double mu_neutral = theta_neut/double(4*N);
  const double mu_del = theta_del/double(4*N);
  const double littler = rho/double(4*N);
  const double m = M/double(4*N);

  gsl_rng * r =  gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(r,seed);

  poptype pop({N}); //Initialize a SINGLE population in the metapop sim!
  std::cerr << typeid(*(pop.gametes.begin())).name() << '\n';
  std::cerr << typeid(pop.gametes).name() << '\n';
  std::function<double(void)> recmap = std::bind(gsl_rng_uniform,r);
  unsigned generation = 0;
  for( unsigned i = 0 ; i < ngens ; ++i,++generation )
    {
      double wbar = KTfwd::sample_diploid(r,
					  &(*pop.gametes.begin()),
					  &pop.diploids[0],
					  &pop.mutations,
					  N,
					  mu_neutral+mu_del,
					  std::bind(KTfwd::infsites(),r,&pop.mut_lookup,
						    mu_neutral,mu_del,[&r](){return gsl_rng_uniform(r);},[&s](){return s;},[&h](){return h;}),
					  std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,
						    &(*pop.gametes.begin()),
						    littler,
						    r,
						    recmap),
					  std::bind(KTfwd::insert_at_end<mtype,mlist>,std::placeholders::_1,std::placeholders::_2),
					  std::bind(KTfwd::insert_at_end<gtype,glist>,std::placeholders::_1,std::placeholders::_2),
					  std::bind(KTfwd::multiplicative_diploid(),std::placeholders::_1,std::placeholders::_2,2.),
					  std::bind(KTfwd::mutation_remover(),std::placeholders::_1,0,2*N));
      KTfwd::remove_fixed_lost(&pop.mutations,&pop.fixations,&pop.fixation_times,&pop.mut_lookup,generation,2*N);
    }

  duplicate_pop(pop.gametes,pop.diploids,pop.mutations);
  std::vector<std::function<double (glist::const_iterator,
				    glist::const_iterator)> > vbf = {
    std::bind(multiplicative_diploid(),std::placeholders::_1,std::placeholders::_2,2.),
    std::bind(multiplicative_diploid_minus(),std::placeholders::_1,std::placeholders::_2,2.)
  };

  unsigned Ns[2] = {N,N};
  double fs[2]={f1,f2};
  for( unsigned i = 0 ; i < ngens2 ; ++i, ++generation )
    {
      std::vector<double> wbars = sample_diploid(r,
						 &pop.gametes,
						 &pop.diploids,
						 &pop.mutations,
						 &Ns[0],
						 mu_neutral + mu_del,
						 std::bind(KTfwd::infsites(),r,&pop.mut_lookup,
							   mu_neutral,mu_del,[&r](){return gsl_rng_uniform(r);},[&s](){return s;},[&h](){return h;}),
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
      remove_fixed_lost(&pop.mutations,&pop.fixations,&pop.fixation_times,&pop.mut_lookup,generation,4*N);
    }
  std::pair< std::vector< std::pair<double,std::string> >,
	     std::vector< std::pair<double,std::string> > > sample1 = ms_sample_separate(r,&pop.diploids[0],n),
    sample2 = ms_sample_separate(r,&pop.diploids[1],n);

  auto neutral = merge(sample1.first,sample2.first,n),
    selected = merge(sample1.second,sample2.second,n);

  //Write the metapop in binary format to outstream
  std::ostringstream outstream;
  KTfwd::write_binary_metapop(&pop.gametes,&pop.mutations,&pop.diploids,
  			      std::bind(KTfwd::mutation_writer(),std::placeholders::_1,std::placeholders::_2),
  			      outstream);

  std::istringstream instream(outstream.str());
  glist_list metapop_gametes2;
  diploid_bucket_vec metapop_diploids2;
  mlist mutations2;

  KTfwd::read_binary_metapop(&metapop_gametes2,&mutations2,&metapop_diploids2,
  			     std::bind(KTfwd::mutation_reader<mtype>(),std::placeholders::_1),
  			     instream);

  //Make sure that the output and the input are the same
  if(pop.mutations.size() != mutations2.size())
    {
      std::cerr << "Error: mutation containers are different sizes. "
  	   << "Line " << __LINE__ << " of " << __FILE__ << '\n';
      std::exit(EXIT_FAILURE);
    }
  auto mitr1 = pop.mutations.begin(),mitr2=mutations2.begin();
  for( ; mitr1 != pop.mutations.end() ; ++mitr1,++mitr2 ) 
    {
      if(mitr1->n != mitr2->n)
  	{
  	  std::cerr << "Error: mutation counts differ. "
  	       << "Line " << __LINE__ << " of " << __FILE__ << '\n';
  	  std::exit(EXIT_FAILURE);
  	}
      if(mitr1->pos != mitr2->pos)
  	{
  	  std::cerr << "Error: mutation positions differ. "
  	       << "Line " << __LINE__ << " of " << __FILE__ << '\n';
  	  std::exit(EXIT_FAILURE);
  	}
      if(mitr1->s != mitr2->s)
  	{
  	  std::cerr << "Error: mutation selection coefficients differ. "
  	       << "Line " << __LINE__ << " of " << __FILE__ << '\n';
  	  std::exit(EXIT_FAILURE);
  	}
      if(mitr1->h != mitr2->h)
  	{
  	  std::cerr << "Error: mutation dominance coefficients differ. "
  	       << "Line " << __LINE__ << " of " << __FILE__ << '\n';
  	  std::exit(EXIT_FAILURE);
  	}
    }

  auto pptr1 = pop.gametes.begin(),pptr2=metapop_gametes2.begin();
  for( unsigned i = 0 ; i < pop.diploids.size() ; ++i,++pptr1,++pptr2 )
    {
      for(unsigned j = 0 ; j < N ; ++j )
  	{
  	  if ( std::distance( pptr1->begin(),pop.diploids[i][j].first ) !=
  	       std::distance( pptr2->begin(),metapop_diploids2[i][j].first ) )
  	    {
  	      std::cerr << "Error: first gametes differ. " 
  			<< "Line " << __LINE__ << " of " << __FILE__ << '\n';
  	      std::exit(EXIT_FAILURE);
  	    }
  	  if ( std::distance( pptr1->begin(),pop.diploids[i][j].second ) !=
  	       std::distance( pptr2->begin(),metapop_diploids2[i][j].second ) )
  	    {
  	      std::cerr << "Error: second gametes differ. " 
  			<< "Line " << __LINE__ << " of " << __FILE__ << '\n';
  	      std::exit(EXIT_FAILURE);
  	    }
  	}
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
