#include <Sequence/SimData.hpp>
#include <Sequence/PolyTableFunctions.hpp>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <utility>
#include <functional>

using namespace std;
using namespace Sequence;

void tpos(SimData & d, const double & x )
{
  transform(d.pbegin(),d.pend(),
	    d.pbegin(),
	    bind2nd(minus<double>(),x));
}

SimData resize( const SimData & d )
{
  SimData rv;

  rv.assign(&*d.pbegin(),d.numsites(),
	    &d[0],2);
  RemoveInvariantColumns(&rv);
  return rv;
}

pair<bool,bool> mono(const SimData & d,
		     const double & lo,
		     const double & hi,
		     const double & L)
{
  unsigned nm1=0,nm2=0;
  for( SimData::const_pos_iterator p = d.pbegin() ; 
       p != d.pend() ; 
       ++p )
    {
      if ( *p*L >= lo - 1. && *p*L < lo )
	{
	  ++nm1;
	}
      else if ( *p*L >= hi - 1. && *p*L < hi )
	{
	  ++nm2;
	}
    }
  return( std::make_pair(!nm1,!nm2) );
}

int main( int argc, char ** argv )
{
  SimData l0,l1,l2;

  while(! cin.eof() )
    {
      cin >> l0 >> l1 >> l2 >> ws;
      if(l0.size()>2)
	{
	  l0=resize(l0);
	}
      if(l1.size()>2)
	{
	  l1=resize(l1);
	}
      if(l1.size()>2)
	{
	  l1=resize(l1);
	}
      tpos(l1,1.);
      tpos(l2,2.);
      pair<bool,bool> m = mono(l0,1.,1000.,1000.);
      cout << m.first << ' ' << m.second << ' ';
      m = mono(l1,1.,1000.,1000.);
      cout << m.first << ' ' << m.second << ' ';
      m = mono(l2,1.,1000.,1000.);
      cout << m.first << ' ' << m.second << '\n';
    }
}
