/*! \include migsel_split_ind_list.cc
  Evolve a constant-N pop, split it off, then evolve both with migration.

  The selection coefficient, s, is treated as -s in deme 2, just for fun.

  Output is a "ms" block for neutral and non-neutral variants.  Ancestral
  population first.
*/
#include <config.h>
#include <numeric>
#include <cmath>
#include <functional>
#include <cassert>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <fwdpp/diploid.hh>

#include <Sequence/SimData.hpp>
#include <Sequence/SimDataIO.hpp> //for writing & reading SimData objects in binary format
#include <Sequence/FST.hpp>
#include <fwdpp/sugar/infsites.hpp>
#include <fwdpp/sugar/infsites.hpp>
#include <fwdpp/sugar/serialization.hpp>
#define METAPOP_SIM
using mtype = KTfwd::mutation;
#include <common_ind.hpp>


using namespace std;
using namespace KTfwd;
using namespace Sequence;

using poptype = metapop_t;
using mcont = poptype::mcont_t;
using gtype = poptype::gamete_t;
using gcont = poptype::gcont_t;

size_t migpop(const size_t & source_pop, gsl_rng * r, const double & mig_prob)
{
  if( gsl_rng_uniform(r) < mig_prob )
    {
      return ! source_pop;
    }
  return source_pop;
}

SimData merge( const std::vector<std::pair<double,std::string> > & sample1,
	       const std::vector<std::pair<double,std::string> > & sample2 ,
	       const unsigned & nsam);

//fitness model details -- s will be treated as -s in population 2
struct multiplicative_diploid_minus
{
  typedef double result_type;
  inline double operator()(const poptype::gamete_t & g1, const poptype::gamete_t & g2,
			   const poptype::mcont_t & mutations,
			   const double scaling = 1.) const
  {
    using mut_t = poptype::mcont_t::value_type;
    return site_dependent_fitness()(g1,g2,mutations,
				    [&](double & fitness,const mut_t & mut)
				    {
				      fitness *= (1. - scaling*mut.s);
				    },
				    [](double & fitness,const mut_t & mut)
				    {
				      fitness *= (1. - mut.h*mut.s);
				    },
				    1.);
  }
};

int main( int argc, char ** argv )
{
  if (argc != 14)
    {
      std::cerr << "Too few arguments.\n"
		<< "Usage: " << argv[0] << " N theta_neutral theta_deleterious rho M s h f1 f2 ngens ngens2 nsam seed\n";
      exit(0);
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
  pop.mutations.reserve(size_t(std::ceil(std::log(2*N)*(theta_neut+theta_del)+0.667*(theta_neut+theta_del))));
  std::function<double(void)> recmap = std::bind(gsl_rng_uniform,r);
  unsigned generation = 0;
  for( unsigned i = 0 ; i < ngens ; ++i,++generation )
    {
      double wbar = KTfwd::sample_diploid(r,
					  pop.gametes,  //non-const reference to gametes
					  pop.diploids[0], //non-const reference to diploids
					  pop.mutations, //non-const reference to mutations
					  pop.mcounts,
					  N,     //current pop size, remains constant
					  mu_neutral+mu_del,    //mutation rate per gamete
					  std::bind(KTfwd::infsites(),std::placeholders::_1,std::placeholders::_2,r,std::ref(pop.mut_lookup),
						    mu_neutral,mu_del,[&r](){return gsl_rng_uniform(r);},[&s](){return s;},[&h](){return h;}),					  
					  //The function to generation recombination positions:
					  std::bind(KTfwd::poisson_xover(),r,littler,0.,1.,
						    std::placeholders::_1,std::placeholders::_2,std::placeholders::_3),
					  std::bind(KTfwd::multiplicative_diploid(),std::placeholders::_1,std::placeholders::_2,
						    std::placeholders::_3,2.),
					  pop.neutral,
					  pop.selected);
      //4*N b/c it needs to be fixed in the metapopulation
      update_mutations(pop.mutations,pop.fixations,pop.fixation_times,pop.mut_lookup,pop.mcounts,generation,4*N);
    }

  //Make an exact copy of the diploids
  pop.diploids.push_back(pop.diploids[0]);

  std::vector<std::function<double(const poptype::gamete_t &,
				   const poptype::gamete_t &,
				   const poptype::mcont_t & )> > vbf;
  vbf.push_back(std::bind(multiplicative_diploid(),std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,2.));
  vbf.push_back(std::bind(multiplicative_diploid_minus(),std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,2.));

  unsigned Ns[2] = {N,N};
  double fs[2]={f1,f2};
  for( unsigned i = 0 ; i < ngens2 ; ++i, ++generation )
    {
      std::vector<double> wbars = sample_diploid(r,
						 pop.gametes,
						 pop.diploids,
						 pop.mutations,
						 pop.mcounts,
						 &Ns[0],
						 mu_neutral + mu_del,
						 std::bind(KTfwd::infsites(),std::placeholders::_1,std::placeholders::_2,r,std::ref(pop.mut_lookup),
							   mu_neutral,mu_del,[&r](){return gsl_rng_uniform(r);},[&s](){return s;},[&h](){return h;}),
						 std::bind(KTfwd::poisson_xover(),r,littler,0.,1.,
							   std::placeholders::_1,std::placeholders::_2,std::placeholders::_3),
						 vbf,
						 std::bind(migpop,std::placeholders::_1,r,m),
						 pop.neutral,
						 pop.selected,
						 &fs[0]);
      //4*N b/c it needs to be fixed in the metapopulation
      update_mutations(pop.mutations,pop.fixations,pop.fixation_times,pop.mut_lookup,pop.mcounts,generation,4*N);
    }
  std::pair< std::vector< std::pair<double,std::string> >,
	     std::vector< std::pair<double,std::string> > > sample1 = ms_sample_separate(r,pop.mutations,pop.gametes,pop.diploids[0],n),
    sample2 = ms_sample_separate(r,pop.mutations,pop.gametes,pop.diploids[1],n);

  auto neutral = merge(sample1.first,sample2.first,n),
    selected = merge(sample1.second,sample2.second,n);

  //Write the metapop in binary format to outstream
  // std::ostringstream outstream;
  // KTfwd::write_binary_metapop(&pop.gametes,&pop.mutations,&pop.diploids,
  // 			      std::bind(KTfwd::mutation_writer(),std::placeholders::_1,std::placeholders::_2),
  // 			      outstream);

  // std::istringstream instream(outstream.str());
  // glist metapop_gametes2;
  // diploid_bucket_vec metapop_diploids2;
  // mlist mutations2;

  // KTfwd::read_binary_metapop(&metapop_gametes2,&mutations2,&metapop_diploids2,
  // 			     std::bind(KTfwd::mutation_reader<mtype>(),std::placeholders::_1),
  // 			     instream);

  // //Make sure that the output and the input are the same
  // if(pop.mutations.size() != mutations2.size())
  //   {
  //     std::cerr << "Error: mutation containers are different sizes. "
  // 	   << "Line " << __LINE__ << " of " << __FILE__ << '\n';
  //     std::exit(EXIT_FAILURE);
  //   }
  // auto mitr1 = pop.mutations.begin(),mitr2=mutations2.begin();
  // for( ; mitr1 != pop.mutations.end() ; ++mitr1,++mitr2 ) 
  //   {
  //     if(mitr1->n != mitr2->n)
  // 	{
  // 	  std::cerr << "Error: mutation counts differ. "
  // 	       << "Line " << __LINE__ << " of " << __FILE__ << '\n';
  // 	  std::exit(EXIT_FAILURE);
  // 	}
  //     if(mitr1->pos != mitr2->pos)
  // 	{
  // 	  std::cerr << "Error: mutation positions differ. "
  // 	       << "Line " << __LINE__ << " of " << __FILE__ << '\n';
  // 	  std::exit(EXIT_FAILURE);
  // 	}
  //     if(mitr1->s != mitr2->s)
  // 	{
  // 	  std::cerr << "Error: mutation selection coefficients differ. "
  // 	       << "Line " << __LINE__ << " of " << __FILE__ << '\n';
  // 	  std::exit(EXIT_FAILURE);
  // 	}
  //     if(mitr1->h != mitr2->h)
  // 	{
  // 	  std::cerr << "Error: mutation dominance coefficients differ. "
  // 	       << "Line " << __LINE__ << " of " << __FILE__ << '\n';
  // 	  std::exit(EXIT_FAILURE);
  // 	}
  //   }

  // auto pptr1 = pop.gametes.begin(),pptr2=metapop_gametes2.begin();
  // for( unsigned i = 0 ; i < pop.diploids.size() ; ++i)//,++pptr1,++pptr2 )
  //   {
  //     for(unsigned j = 0 ; j < N ; ++j )
  // 	{
  // 	  if ( std::distance( pptr1,pop.diploids[i][j].first ) !=
  // 	       std::distance( pptr2,metapop_diploids2[i][j].first ) )
  // 	    {
  // 	      std::cerr << "Error: first gametes differ. " 
  // 			<< "Line " << __LINE__ << " of " << __FILE__ << '\n';
  // 	      std::exit(EXIT_FAILURE);
  // 	    }
  // 	  if ( std::distance( pptr1,pop.diploids[i][j].second ) !=
  // 	       std::distance( pptr2,metapop_diploids2[i][j].second ) )
  // 	    {
  // 	      std::cerr << "Error: second gametes differ. " 
  // 			<< "Line " << __LINE__ << " of " << __FILE__ << '\n';
  // 	      std::exit(EXIT_FAILURE);
  // 	    }
  // 	}
  //   }
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
