#include <fwdpp/diploid.hh>
#include <Sequence/SimData.hpp>
#include <boost/unordered_set.hpp>

#include <boost/container/list.hpp>
#include <boost/container/vector.hpp>
#include <boost/pool/pool_alloc.hpp>
#include <boost/function.hpp>
#include <vector>
#include <list>

using namespace std;
using namespace Sequence;
using namespace KTfwd;

struct mutation_with_age : public KTfwd::mutation_base
{
  unsigned g;
  double s,h;
  mutation_with_age(const unsigned & __o,const double & position, const unsigned & count, const bool & isneutral = true)
    : KTfwd::mutation_base(position,count,isneutral),g(__o),s(0.),h(0.)
  {	
  }
};

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


typedef boost::unordered_set<double,boost::hash<double>,KTfwd::equal_eps > lookup_table_type;
mutation_with_age neutral_mutations_inf_sites(gsl_rng * r,const unsigned * generation,mlist * mutations,
					      lookup_table_type * lookup, const double & beg)
{
  double pos = gsl_ran_flat(r,beg,beg+1.);
  while( lookup->find(pos) != lookup->end() ) 
    { 
      pos = gsl_ran_flat(r,beg,beg+1.);
    }
  lookup->insert(pos);
  assert(std::find_if(mutations->begin(),mutations->end(),boost::bind(KTfwd::mutation_at_pos(),_1,pos)) == mutations->end());
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
  const unsigned seed = (argc==1) ? 0 : atoi(argv[argument++]);        //Random number seed
  const unsigned N = 1000; //Population size
  const double L = 1000; //Locus length in psuedo-sites a-la ms.
  const double theta = 1;
  const double rho = 1;
  const unsigned ngens = 10*N;
  const unsigned samplesize = 2;
  int nreps = 1000;

  const double mu = theta/double(4*N);                 //per-gamete mutation rate
  const double littler = rho/double(4*N);
  
  //Initiate random number generation system
  gsl_rng * r =  gsl_rng_alloc(gsl_rng_ranlxs2);
  gsl_rng_set(r,seed);

  unsigned twoN = 2*N;

  //recombination map is uniform[0,1) for all three loci
  boost::function<double(void)> recmap = boost::bind(gsl_rng_uniform,r),
    recmap2 = boost::bind(gsl_ran_flat,r,1.,2.),
    recmap3 = boost::bind(gsl_ran_flat,r,2.,3.);

  //need vectors of recombination maps, mutation policies, and fitness models
  std::vector< boost::function<unsigned(glist::iterator &, glist::iterator &)> > recpols(3);
  std::vector< boost::function<mtype(mlist *)> > mmodels(3);

  const double rbw[2] = {0.001,0.0001};
  while(nreps--)
    {
      std::vector< glist > gametes (3, glist(1,gtype(twoN) ));
      std::vector< std::pair< glist::iterator, glist::iterator > > idip(3);
      idip[0] = std::make_pair( gametes[0].begin(),gametes[0].begin() );
      idip[1] = std::make_pair( gametes[1].begin(),gametes[1].begin() );
      idip[2] = std::make_pair( gametes[2].begin(),gametes[2].begin() );

      std::vector< std::vector< std::pair< glist::iterator, glist::iterator > > > diploids(N,idip);

      mlist mutations;  
      std::vector<mtype> fixations;  
      std::vector<unsigned> fixation_times;
      unsigned generation=0;
      double wbar;
      lookup_table_type lookup;

      recpols[0] = boost::bind(KTfwd::genetics101(),_1,_2,&gametes[0],littler,r,recmap);
      recpols[1] = boost::bind(KTfwd::genetics101(),_1,_2,&gametes[1],littler,r,recmap2);
      recpols[2] = boost::bind(KTfwd::genetics101(),_1,_2,&gametes[2],littler,r,recmap3);
      mmodels[0] = boost::bind(neutral_mutations_inf_sites,r,&generation,_1,&lookup,0.);
      mmodels[1] = boost::bind(neutral_mutations_inf_sites,r,&generation,_1,&lookup,1.);
      mmodels[2] = boost::bind(neutral_mutations_inf_sites,r,&generation,_1,&lookup,2.);
      for( generation = 0; generation < ngens; ++generation )
      	{
	  KTfwd::sample_diploid( r,
	  			 &gametes,
	  			 &diploids,
	  			 &mutations,
	  			 N,
	  			 N,
	  			 mu,
	  			 mmodels,
	  			 recpols,
	  			 rbw,
	  			 boost::bind(KTfwd::insert_at_end<mtype,mlist>,_1,_2),
	  			 boost::bind(KTfwd::insert_at_end<gtype,glist>,_1,_2),
	  			 boost::bind(no_selection_multi(),_1),
	  			 boost::bind(KTfwd::mutation_remover(),_1,0,2*N),
	  			 0.);
      	  KTfwd::remove_fixed_lost(&mutations,&fixations,&fixation_times,&lookup,generation,2*N);
	}
      vector< vector< std::pair<double,string> > > gam_sample = ms_sample(r,&diploids,samplesize,true);
									  
    }
  return 0;
}
