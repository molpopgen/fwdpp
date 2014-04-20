#include <Sequence/SimData.hpp>
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

pair<bool,bool> mono(const SimData & d,
		     const double & lo,
		     const double & hi)
{
  unsigned c1 = count_if(d.pbegin(),d.pend(),bind2nd(less<double>(),lo));
  unsigned c2 = count_if(d.pbegin(),d.pend(),bind2nd(greater_equal<double>(),hi));
  return make_pair(c1,c2);
}
int main( int argc, char ** argv )
{
  SimData l0,l1,l2;

  while(! cin.eof() )
    {
      cin >> l0 >> l1 >> l2 >> ws;
      tpos(l1,1.);
      tpos(l2,2.);
      pair<bool,bool> m = mono(l0,1./1000.,1.-1./1000.);
      cout << m.first << ' ' << m.second << ' ';
      m = mono(l1,1./1000.,1.-1./1000.);
      cout << m.first << ' ' << m.second << ' ';
      m = mono(l2,1./1000.,1.-1./1000.);
      cout << m.first << ' ' << m.second << '\n';
    }
}
