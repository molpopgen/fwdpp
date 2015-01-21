#include <Sequence/SimData.hpp>
#include <iostream>

using namespace std;
using namespace Sequence;

int main(int argc, char **argv)
{
  int rv;
  unsigned nmono=0,nruns=0;
  SimData d;
  
  const unsigned site1 = atoi(argv[1]);
  const unsigned site2 = atoi(argv[2]);
  const unsigned nsites = atoi(argv[3]);

  unsigned low=0,ttl=0;
  while( (rv=d.fromfile(stdin)) != EOF )
    {
      ++nruns;
      unsigned nm1=0,nm2=0;
      ttl += d.numsites();
      for( unsigned i=0; i < d.numsites() ; ++i )
        {
          if(d.position(i) <= 0.5)
            ++low;
        }
      if(d.numsites()==0) ++nmono;
      else
        {
          for( unsigned i=0; i < d.numsites() ; ++i )
            {
              if( d.position(i)*double(nsites) >= double(site1-1) && 
                  d.position(i)*double(nsites) < double(site1) )
                {
                  ++nm1;
                }
              else if( d.position(i)*double(nsites) >= double(site2-1) && 
                       d.position(i)*double(nsites) < double(site2) )
                {
                  ++nm2;
                }
            }
          if(nm1==0 && nm2 == 0)
            ++nmono;
        }
    }
  cout << nmono << ' ' << nruns << ' ' << double(nmono)/double(nruns) << ' ' << low << ' ' << ttl << '\n';
}
