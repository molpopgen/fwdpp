/*! \include diploid_ind_2locus.cc
  Simple example of a two-locus simulation using the multilocus API in fwdpp.
*/

#include <fwdpp/diploid.hh>
#include <Sequence/SimData.hpp>
#include <vector>
#include <list>
#include <sstream>

/*
  We define a new mutation type, derived from the base class KTfwd::mutation_base.
  This type adds a selection coefficient (s), dominance coefficient (h),
  and the generation when the mutation first appeared in the population (g)
*/
struct mutation_with_age : public KTfwd::mutation_base
{
  unsigned g;
  double s,h;
  /*
    The constructor initializes all class data, including that of the base class via a constructor
    call to the base class.
  */
  mutation_with_age(const unsigned & __o,const double & position, const unsigned & count, const bool & isneutral = true)
    : KTfwd::mutation_base(position,count,isneutral),g(__o),s(0.),h(0.)
  {	
  }
};

using mtype = mutation_with_age;
//function object to write mutation data in binary format
struct mwriter
{
  typedef void result_type;
  result_type operator()( const mutation_with_age & m, std::ostringstream & buffer ) const
  {
    unsigned u = m.n;
    buffer.write( reinterpret_cast< char * >(&u),sizeof(unsigned) );
    u = m.g;
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
  typedef mutation_with_age result_type;
  result_type operator()( std::istream & in ) const
  {
    unsigned n;
    in.read( reinterpret_cast< char * >(&n),sizeof(unsigned) );
    unsigned g;
    in.read( reinterpret_cast< char * >(&g),sizeof(unsigned) );
    bool neut;
    in.read( reinterpret_cast< char * >(&neut),sizeof(bool) );
    double pos;
    in.read( reinterpret_cast< char * >(&pos),sizeof(double) );
    double s;
    in.read( reinterpret_cast< char * >(&s),sizeof(double) );
    double h;
    in.read( reinterpret_cast< char * >(&h),sizeof(double) );
    return result_type(g,pos,n,neut);
  }
};
#include <common_ind.hpp>

/*
  This function generates neutral mutations under the infinitely-many sites assumption
*/
mutation_with_age neutral_mutations_inf_sites(gsl_rng * r,const unsigned * generation,mlist * mutations,
					      lookup_table_type * lookup, const double & beg)
{
  //Generate new mutation position on the interval [0,1)
  double pos = gsl_ran_flat(r,beg,beg+1.);
  while( lookup->find(pos) != lookup->end() ) //make sure it doesn't exist in the population
    { 
      pos = gsl_ran_flat(r,beg,beg+1.);
    }
  //update the lookup table
  lookup->insert(pos);

  //In absence of DEBUG, make sure lookup table is working
  assert(std::find_if(mutations->begin(),mutations->end(),std::bind(KTfwd::mutation_at_pos(),std::placeholders::_1,pos)) == mutations->end());

  //return constructor call to mutation type
  return mutation_with_age(*generation,pos,1,true);
}
 
//Fitness function
struct no_selection_multi
{
  typedef double result_type;
  template< typename iterator_type >
  inline double operator()( const std::vector< std::pair<iterator_type,iterator_type> > & diploid ) const
  {
    return 1.;
  }
};

int main(int argc, char ** argv)
{
  int argument=1;
  if (argc != 8 )
    {
      std::cerr << "Incorrect number of arguments.\n"
		<< "Usage:\n"
		<< argv[0] << " N theta rho ngens n nreps seed\n"
		<< "Where:\n"
		<< "N = population size (number of diploids)\n"
		<< "theta = 4Nu, the scaled neutral mutation rate\n"
		<< "rho = 4Nr, the scale recombination rate\n"
		<< "ngens = the number of generations to simulate\n"
		<< "n = the sample size to pull from the population at the end of each simulated replicate\n"
		<< "nreps = the number of replicates to simulated\n"
		<< "seed = seed value for random number generations\n";
	std::exit(0);
    }
  const unsigned N = atoi(argv[argument++]);           //Number of diploids
  const double theta = atof(argv[argument++]);         //4*n*mutation rate.  Note: mutation rate is per REGION, not SITE!!
  const double rho = atof(argv[argument++]);           //4*n*recombination rate.  Note: recombination rate is per REGION, not SITE!!
  const unsigned ngens = atoi(argv[argument++]);       //Number of generations to simulate
  const unsigned samplesize1 = atoi(argv[argument++]); //Sample size to draw from the population
  int nreps = atoi(argv[argument++]);                  //Number of replicates to simulate
  const unsigned seed = atoi(argv[argument++]);        //Random number seed

  const std::vector<double> mu(2,theta/double(4*N));   //per-gamete mutation rate per locus
  
  /*
    littler r is the recombination rate per region per generation.

    For individual simulation (UNLIKE GAMETE-BASED SIMS!!!),
    r = rho/(4N)
  */
  const double littler = rho/double(4*N);
  
  //Initiate random number generation system
  gsl_rng * r =  gsl_rng_alloc(gsl_rng_ranlxs2);
  gsl_rng_set(r,seed);

  unsigned twoN = 2*N;

  std::function<double(void)> recmap = std::bind(gsl_rng_uniform,r),
    recmap2 = std::bind(gsl_ran_flat,r,1.,2.);

  const double rbw = 0.;//0.1;
  while(nreps--)
    {
      //the population begins with 1 gamete with no mutations amd initial count 2N
      std::vector< glist > gametes (2, glist(1,gtype(twoN) ));
      //Establish the list of diploids.  Here, there is only 1 gamete, so it is easy:
      std::vector< std::pair< glist::iterator, glist::iterator > > idip(2);
      idip[0] = std::make_pair( gametes[0].begin(),gametes[0].begin() );
      idip[1] = std::make_pair( gametes[1].begin(),gametes[1].begin() );
      std::vector< std::vector< std::pair< glist::iterator, glist::iterator > > > diploids(N,idip);

      mlist mutations;  //the population is devoid of mutations
      std::vector<mtype> fixations;  //store mutation that fix.  Passed to KTfwd::remove_fixed_lost
      std::vector<unsigned> fixation_times; //store when they fix.  Passed to KTfwd::remove_fixed_lost
      unsigned generation=0;
      double wbar;
      lookup_table_type lookup;  //this is our lookup table for the mutation model
      //within-locus recombination policies -- one per locus
      auto recpol0 = std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,&gametes[0],littler,r,recmap);
      auto recpol1 = std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,&gametes[1],littler,r,recmap2);
      std::vector< decltype(recpol0) > recpols{ recpol0 , recpol1 };

      //mutation policies -- 1 per locus
      auto mmodel0 = std::bind(neutral_mutations_inf_sites,r,&generation,std::placeholders::_1,&lookup,0.);
      auto mmodel1 = std::bind(neutral_mutations_inf_sites,r,&generation,std::placeholders::_1,&lookup,1.);
      std::vector< decltype(mmodel0) > mmodels { mmodel0, mmodel1 };
      for( generation = 0; generation < ngens; ++generation )
      	{
      	  //Iterate the population through 1 generation
	  KTfwd::sample_diploid( r,
	  			 &gametes,
	  			 &diploids,
	  			 &mutations,
	  			 N,
	  			 N,
				 &mu[0],
	  			 mmodels,
	  			 recpols,
	  			 &rbw,
				 [](gsl_rng * __r, const double __d){ return gsl_ran_poisson(__r,__d); },
	  			 std::bind(KTfwd::insert_at_end<mtype,mlist>,std::placeholders::_1,std::placeholders::_2),
	  			 std::bind(KTfwd::insert_at_end<gtype,glist>,std::placeholders::_1,std::placeholders::_2),
	  			 std::bind(no_selection_multi(),std::placeholders::_1),
	  			 std::bind(KTfwd::mutation_remover(),std::placeholders::_1,0,2*N),
	  			 0.);
	  assert( check_sum(gametes[0],twoN) );
	  assert( check_sum(gametes[1],twoN) );
      	  KTfwd::remove_fixed_lost(&mutations,&fixations,&fixation_times,&lookup,generation,2*N);
	}
      unsigned nm1=0,nm2=0;
      for(mlist::const_iterator i = mutations.begin() ; i != mutations.end() ; ++i )
	{
	  if ( i->pos < 1.)++nm1;
	  else ++nm2;
	}
      std::cout << nm1 << '\t' << nm2 << '\n';
    }
  return 0;
}
