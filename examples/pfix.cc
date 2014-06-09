/*
  Estimate fixation probability under single-site model with additivity.  P(fix) should be \approx 2s for small s.
*/

#include <fwdpp/diploid.hh>
#include <iostream>
#include <cstdlib>
#include <boost/bind.hpp>
typedef KTfwd::mutation mtype;
typedef std::list<mtype> mlist;
typedef KTfwd::gamete_base<mtype,mlist> gtype;
typedef std::list<gtype> glist;
typedef std::vector< std::pair<glist::iterator,
			       glist::iterator> > dipvec;
using namespace std;
using namespace KTfwd;

struct no_mutation
{
  typedef mtype result_type;
  template<typename T>
  inline mtype operator()(T & g)const
  {
    abort();//must never be called!
    return mtype(0,0,1,1);
  }
};

struct no_recombination
{
  inline double operator()(void) const
  {
    return 0.;
  }
};
int main( int argc, char ** argv )
{
  int argument=1;
  const unsigned N = atoi(argv[argument++]);
  const double s = atof(argv[argument++]);
  const unsigned nreps = atoi(argv[argument++]);
  const unsigned seed = atoi(argv[argument++]);

  gsl_rng * r =  gsl_rng_alloc(gsl_rng_ranlxs2);
  gsl_rng_set(r,seed);

  unsigned fixations = 0,rep=0;
  while( rep++ < nreps )
    {
      mlist mutations;
      glist gametes(1,gtype(2*N));

      dipvec diploids(2*N,make_pair(gametes.begin(),gametes.begin()));

      //manually add a selected mutation to population at position 0.5 and dominance of 1
      mlist::iterator mitr = mutations.insert(mutations.end(),
					      mutation(0.5,s,1,1.));
      //make a new gamete containing this new mutation
      gtype mutant(1);
      mutant.smutations.push_back(mitr);
      //make first gamete of 1st diploid the mutant
      diploids[0].first = gametes.insert(gametes.end(),mutant);

      while(mitr->n > 0 && mitr->n < 2*N)
	{
	  sample_diploid(r,
			 &gametes,
			 &diploids,
			 &mutations,
			 N,
			 0.,//mutation rate
			 boost::bind( no_mutation(),_1 ),
			 boost::bind(KTfwd::genetics101(),_1,_2,
				     &gametes,
				     0.,
				     r,
				     no_recombination()),
			 boost::bind(KTfwd::insert_at_end<mtype,mlist>,_1,_2),
			 boost::bind(KTfwd::insert_at_end<gtype,glist>,_1,_2),
			 boost::bind(KTfwd::multiplicative_diploid(),_1,_2,2.),
			 boost::bind(KTfwd::mutation_remover(),_1,0,2*N));
	}
      fixations += (mitr->n == 2*N) ? 1 : 0;
    }
  cout << double(fixations)/double(nreps) << '\n';
}
