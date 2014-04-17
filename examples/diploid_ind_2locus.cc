/*! \include diploid_ind.cc
  Simulate a single, finite Wright-Fisher diploid population with mutation, recombination, and no selection.

  This is an individual-based implementation of the gamete-based program diploid.cc

  This program illustrates many features of fwdpp:
  1.  Custom mutation classes
  2.  Implementing a mutation model (infinitely-many sites)
  3.  Iterating a population through its life cycle
  4.  Outputting a sample in "ms" format
*/

#include <fwdpp/diploid.hh>
#include <Sequence/SimData.hpp>
#include <boost/unordered_set.hpp>

#include <boost/container/list.hpp>
#include <boost/container/vector.hpp>
#include <boost/pool/pool_alloc.hpp>
#include <boost/function.hpp>
#include <vector>
#include <list>
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

//Simplify our lives a bit with typedefs

/*
  NOTE: we are using boost containers with boost memory allocators.  
  My experience on my systems is that each of those objects is worth a 5-10% 
  run time speedup compare to the same objects from namespace std.

  fwdpp is compatible with ANY container system as long as the
  one containing gametes conforms to the behavior of std::vector and 
  the one containing mutations conforms to the behavior of std::list 
  (especially w/regards to no pointer invalidation upon insertion/delete!!)
*/

//compiling the code with -DUSE_STANDARD_CONTAINERS will use std::vector and std::list instead of the boost alternatives
typedef mutation_with_age mtype;
#ifndef USE_STANDARD_CONTAINERS
typedef boost::pool_allocator<mtype> mut_allocator;
typedef boost::container::list<mtype,mut_allocator > mlist;
typedef KTfwd::gamete_base<mtype,mlist> gtype;
typedef boost::pool_allocator<gtype> gam_allocator;
typedef boost::container::list<gtype,gam_allocator > glist;
#else
typedef std::list<mtype> mlist;
typedef KTfwd::gamete_base<mtype,mlist> gtype;
typedef std::list<gtype > list;
#endif


/*
  We wish to do mutations under the infinitely-many sites assumption.  That means that
  a new mutation cannot appear at any previously-mutated site.  Here, we cheat a little
  and do not allow mutations at any sites that are currently polymorphic.

  We accomplish this via a lookup table of the current mutation positions.  The function object
  KTfwd::equal_eps is used as a replacement for std::operator==(double,double) in order to ensure
  that values differing by <= DBL_EPSILON (~10^-17 on most systems) are not allowed, as they cause problems.
*/
typedef boost::unordered_set<double,boost::hash<double>,KTfwd::equal_eps > lookup_table_type;

/*
  This function generates neutral mutations under the infinitely-many sites assumption
*/
mutation_with_age neutral_mutations_inf_sites(gsl_rng * r,const unsigned * generation,mlist * mutations,
					      lookup_table_type * lookup, const double & beg)
{
  //Generate new mutation position on the interval [0,1)
  //double pos = gsl_rng_uniform(r);
  double pos = gsl_ran_flat(r,beg,beg+1.);
  /*
    An alternative implementation of the while loop below would be:
    while( std::find_if( mutations->begin(),mutations->end(),boost::bind(KTfwd::mutation_at_pos(),_1,pos) != mutations->end()) )
    {
    pos = gsl_rng_uniform(r);
    }

    However, that operation is typically much slower, esp. as the mutation rate gets higher
  */
  while( lookup->find(pos) != lookup->end() ) //make sure it doesn't exist in the population
    { 
      pos = gsl_ran_flat(r,beg,beg+1.);
      //pos = gsl_rng_uniform(r);  //if it does, generate a new one
    }
  //update the lookup table
  lookup->insert(pos);

  /*
    The program will abort if the following conditions are true:
    1.  The mutation position that we generated does indeed exist in the population (implying that our hashing scheme has failed).
    2.  AND, the program was NOT compiled with -DNDEBUG

    fwdpp makes extensive internal use of the C-language assert macro (#include <cassert>).  The macro is enabled if a program is 
    compiled like this:
    c++ -o program program.cc

    It is disabled if a program is compiled like this:
    c++ -DNDEBUG program program.cc

    Thus, a programmer may make use of the preprocessor macro NDEBUG to write detailed debugging routines that are easily disabled.
    
    For example:
    #ifndef NDEBUG
    bool condition = false;
    //do some complex checking of your data here.  If all is okay, set condition = true;
    #endif
    assert( condition == true );

    That check may be expensive at run-time, but is easily disabled by recompiling using -DNDEBUG

    An example of such a debugging routine is KTfwd::check_sum
  */
  assert(std::find_if(mutations->begin(),mutations->end(),boost::bind(KTfwd::mutation_at_pos(),_1,pos)) == mutations->end());

  //return constructor call to mutation type
  return mutation_with_age(*generation,pos,1,true);
}
 
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
  const unsigned N = atoi(argv[argument++]);           //Number of diploids
  const double theta = atof(argv[argument++]);         //4*n*mutation rate.  Note: mutation rate is per REGION, not SITE!!
  const double rho = atof(argv[argument++]);           //4*n*recombination rate.  Note: recombination rate is per REGION, not SITE!!
  const unsigned ngens = atoi(argv[argument++]);       //Number of generations to simulate
  const unsigned samplesize1 = atoi(argv[argument++]); //Sample size to draw from the population
  int nreps = atoi(argv[argument++]);                  //Number of replicates to simulate
  const unsigned seed = atoi(argv[argument++]);        //Random number seed

  //const double mu = theta/double(4*N);                 //per-gamete mutation rate
  const std::vector<double> mu(2,theta/double(4*N));
  
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

  //recombination map is uniform[0,1)
  boost::function<double(void)> recmap = boost::bind(gsl_rng_uniform,r),
    recmap2 = boost::bind(gsl_ran_flat,r,1.,2.);

  //need vectors of recombination maps, mutation policies, and fitness models
  std::vector< boost::function<unsigned(glist::iterator &, glist::iterator &)> > recpols(2);
  std::vector< boost::function<mtype(mlist *)> > mmodels(2);

  const double rbw = 0.;//0.1;
  while(nreps--)
    {
      //the population begins with 1 gamete with no mutations amd initial count 2N
      //glist gametes(1,gtype(twoN));
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
      //rec policies
      recpols[0] = boost::bind(KTfwd::genetics101(),_1,_2,&gametes[0],littler,r,recmap);
      recpols[1] = boost::bind(KTfwd::genetics101(),_1,_2,&gametes[1],littler,r,recmap2);
      mmodels[0] = boost::bind(neutral_mutations_inf_sites,r,&generation,_1,&lookup,0.);
      mmodels[1] = boost::bind(neutral_mutations_inf_sites,r,&generation,_1,&lookup,1.);
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
	  			 //mu,
	  			 mmodels,
	  			 recpols,
	  			 &rbw,
	  			 boost::bind(KTfwd::insert_at_end<mtype,mlist>,_1,_2),
	  			 boost::bind(KTfwd::insert_at_end<gtype,glist>,_1,_2),
	  			 boost::bind(no_selection_multi(),_1),
	  			 boost::bind(KTfwd::mutation_remover(),_1,0,2*N),
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
