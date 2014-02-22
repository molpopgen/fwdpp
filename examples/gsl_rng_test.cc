/*
  How often does gsl_rng_uniform return the same value?

  Usage: ./gsl_rng_test seed nreps
*/

#include <iostream>
#include <cstdlib>

#include <boost/unordered_set.hpp>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;

typedef boost::unordered_set<double,boost::hash<double> > lookup_table_type;

int main( int argc, char ** argv )
{
  unsigned seed = atoi(argv[1]);
  unsigned reps = atoi(argv[2]);

  gsl_rng * r =  gsl_rng_alloc(gsl_rng_ranlxs2);
  gsl_rng_set(r,seed);

  double sum = 0.;
  for( unsigned i = 0 ; i < reps ; ++i )
    {
      lookup_table_type l;
      long unsigned NCALLS = 0;
      bool done = false;
      while(! done )
	{
	  double p = gsl_rng_uniform(r);
	  ++NCALLS;
	  if( l.find(p) != l.end() )
	    {
	      done = true;
	      break;
	    }
	  else
	    {
	      l.insert(p);
	    }
	}
      sum += double(NCALLS);
    }
  cout << "On average, you get the same value every " << sum/double(reps) << " calls to gsl_rng_uniform\n";
}
