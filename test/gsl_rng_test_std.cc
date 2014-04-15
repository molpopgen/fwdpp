/*
  How often does gsl_rng_uniform return the same value?

  Usage: ./gsl_rng_test seed nreps
*/

#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <vector>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;

int main( int argc, char ** argv )
{
  unsigned seed = atoi(argv[1]);
  unsigned reps = atoi(argv[2]);

  gsl_rng * r =  gsl_rng_alloc(gsl_rng_ranlxs2);
  gsl_rng_set(r,seed);

  double sum = 0.;
  for( unsigned i = 0 ; i < reps ; ++i )
    {
      vector<double> l;
      long unsigned NCALLS = 0;
      bool done = false;
      while(! done )
	{
	  double p = gsl_rng_uniform(r);
	  ++NCALLS;
	  if( find( l.begin(),l.end(),p ) != l.end() )
	    {
	      done = true;
	      break;
	    }
	  else
	    {
	      l.push_back(p);
	    }
	}
      sum += double(NCALLS);
    }
  cout << "On average, you get the same value every " << sum/double(reps) << " calls to gsl_rng_uniform\n";
}
