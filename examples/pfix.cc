/*
  Estimate fixation probability under single-site model with additivity.  P(fix) should be \approx 2s for small s.

  This simulation is equvalent to the following R code:

  pfixWF=function(N,s,nreps)
  {
    nfixed=0
    for( i in 1:nreps )
      {
        q=1/(2*N)
        while(q>0&&q<1)
        {
          p = 1-q
          wbar = p^2 + 2*p*q*(1+s) + (1+2*s)*q^2
          Eqprime = q + p*q*s/wbar
          q = rbinom(1,2*N,Eqprime)/(2*N)
        }
        if(q==1)
          {
            nfixed=nfixed+1
          }
      }
    return(nfixed/nreps)
  }

  However, this simulation will be SLOWER than the R code because it is doing a lot more each generation.

  In essence, using fwdpp like this is like diggin a fence post hole with dynamite.
*/

#include <fwdpp/diploid.hh>
#include <iostream>
#include <cstdlib>
#include <boost/bind.hpp>

/*
  typedefs simplify life.

  These use containers from namespace std.

  Containers from boost are often faster 
  and use the same syntax, but we'll keep
  things "std" in this example.
 */

/*
  Use the basic mutation type from fwdpp,
  including s and h as selection coefficient
  and dominance
*/
typedef KTfwd::mutation mtype;
/*
  Store mutations in doubly-linked list
*/
typedef std::list<mtype> mlist;
/*
  Gametes container iterators derived from mlist,
  so they need to know our mutation type.
 */
typedef KTfwd::gamete_base<mtype,mlist> gtype;
/*
  glist is our doubly-linked list of gametes
 */
typedef std::list<gtype> glist;
/*
  Finally, dipvec is a vector of pairs of iterators
  derived from our glist
 */
typedef std::vector< std::pair<glist::iterator,
			       glist::iterator> > dipvec;

//The relevant namespaces
using namespace std;
using namespace KTfwd;

/*
  This example has a mutation rate of 0. 

  However, a mutation policy MUST
  return a mutation object.

  This policy's operator() returns
  a mutation at position 0, s = 0,
  count = 1, and dominance = 1.
 */
struct no_mutation
{
  typedef mtype result_type;
  //A pointer to a list of mutations MUST be passed to the mutation model
  inline mtype operator()(mlist * mutations)const 
  {
    abort();//must never be called!
    return mtype(0,0,1,1);
  }
};

/*
  This example has no recombination.

  However, fwdpp requires a genetic map
  policy that returns a double, representing
  the position of an xover event.

  This is that policy for this sim.
 */
struct no_recombination
{
  inline double operator()(void) const
  {
    abort();//this function must never be called!
    return 0.;
  }
};

int main( int argc, char ** argv )
{
  if ( argc != 5)
    {
      cerr << "Too few arguments.\n"
	   << "Usage: pfix N s nreps seed\n";
      exit(10);
    }
  int argument=1;
  const unsigned N = atoi(argv[argument++]);
  const double s = atof(argv[argument++]);
  const unsigned nreps = atoi(argv[argument++]);
  const unsigned seed = atoi(argv[argument++]);

  /*
    fwdpp uses GSL's random number generation scheme.
    See the following for documentation:
    http://www.gnu.org/software/gsl/manual/html_node/Random-Number-Generation.html#Random-Number-Generation
    http://www.gnu.org/software/gsl/manual/html_node/Random-Number-Distributions.html#Random-Number-Distributions

    Requires -lgsl -lgslcblas at link time.

    Note: In C, an additional -lm would be required AFTER linking the GSL librarys.  However, the -lm
    is implicit in C++ because the standard C++ library requires the math library.
   */
  gsl_rng * r =  gsl_rng_alloc(gsl_rng_ranlxs2);
  gsl_rng_set(r,seed);

  //Now, we cound the number of fixations in the next nreps simulations
  unsigned fixations = 0,rep=0;
  while( rep++ < nreps )
    {
      //Create an empty mutation list
      mlist mutations;
      //The population begins with 1 gamete and 2N copies of that gamete
      glist gametes(1,gtype(2*N));

      //The population begins with N diploids.  Each diploid contains two copies of the initial gamete
      dipvec diploids(N,make_pair(gametes.begin(),gametes.begin()));

      /*
	The next several operations give some insight into what the internals of the library are doing.
	Typically, we'd pick a gamete, decide how many mutations it gets, then apply mutation model that
	many times to the gamete, etc.  Instead, we will manually add our selected mutation to a 
	randomly-chosen gamete.
      */

      /*
	Manually add a selected mutation to population at position 0.5 and dominance of 1.
       The return value of an insert into a linked list is an iterator recording the position
       in the list of where the insertion took place.  We need that info.
      */
      mlist::iterator mitr = mutations.insert(mutations.end(),
					      mutation(0.5,s,1,1.));
      //make a new gamete containing this new mutation
      gtype mutant(1);

      //The mutation is selected, so we add it to the selected mutations container of our new gamete
      mutant.smutations.push_back(mitr);


      //make a random gamete the mutant.  Pick individual first, then gamete
      unsigned ind = unsigned(gsl_ran_flat(r,0.,double(N)));
      unsigned gam = (gsl_rng_uniform(r)<=0.5) ? 0 : 1;
      if(!gam)
	{
	  diploids[ind].first = gametes.insert(gametes.end(),mutant);
	}
      else
	{
	  diploids[ind].second = gametes.insert(gametes.end(),mutant);
	}

      //We're now beyond the part where we are doing "internal" operations manually.

      while(mitr->n > 0 && mitr->n < 2*N)//While our mutation is still segregating...
	{
	  //...we sample the next generation.

	  /*
	    sample_diploid is the main function for individual-based simulations.
	  */
	  sample_diploid(r,
			 &gametes,
			 &diploids,
			 &mutations,
			 N,                                                         //If we wished pop. size to change, we'd pass N,N2 here.  For models where N in next generation changes over time, N2 would be updated accordingly each generation at the end of this block.
			 0.,                                                        //mutation rate = 0.
			 boost::bind( no_mutation(),_1 ),                           //Our mutation model.  The _1 is a placeholder for a pointer to the list of mutations.
			 boost::bind(KTfwd::genetics101(),_1,_2,                    //genetics101 is the basic model of crossing over in fwdpp.  If you want different models, write a replacement policy.  _1 and _2 are place holders for the two gametes that will recombine.
				     &gametes,
				     0.,                                            //The rec. rate is 0
				     r,
				     no_recombination()),                           //This is our genetic map
			 boost::bind(KTfwd::insert_at_end<mtype,mlist>,_1,_2),      //This policy defines how new mutations are added to the mlist.  They are inserted at the end.  This will never be called (as mut. rate = 0), but it is required.
			 boost::bind(KTfwd::insert_at_end<gtype,glist>,_1,_2),      //This policy defines how new gametes are added to the glist.  They are inserted at the end.  This will never be called (no mutation nor recombination), but it is required.
			 boost::bind(KTfwd::multiplicative_diploid(),_1,_2,2.),     //The fitness model.  Fitnesses are 1, 1+hs, 1+scaling*s per site, and multiplicative over sites.  Scaling = 2 and h = 1 (as used here) gives 1, 1+s, 1+2s, which is genic selection.
			 boost::bind(KTfwd::mutation_remover(),_1,0,2*N));          //A mutation that hits a count of either 0 or 2N is removed from a gamete.

	  /*
	    Typically, for a "pop-gen" simulation where relative fitness is all that matters,
	    after sampling the daughter generation, we would want to remove 
	    fixed/lost mutations from the mlist via a call to remove_fixed_lost.

	    However, this is one of those simultations where we don't want to do that
	    b/c we want to know the final fate fo the mutation.

	    Another example of where we would not want to remove fixed mutations would
	    be for quantitative trait sims and/or sims involving evolution to new phenotypic
	    optimum. In those cases, fixed mutations affect trait values and should be kept around.
	   */
	}
      //increment our fixations counter
      fixations += (mitr->n == 2*N) ? 1 : 0;
    }
  //Write fixation probability to STDOUT
  cout << double(fixations)/double(nreps) << '\n';
}
