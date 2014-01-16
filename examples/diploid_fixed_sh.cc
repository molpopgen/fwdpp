/*
  \include diploid_fixed_sh.cc
  
  Same as diploid.cc, but extended to add a mutation rate to beneficial
  alleles with selection coefficient s and dominance h, both provided on the command line
 */

#include <fwdpp/diploid.hh>
#include <boost/unordered_set.hpp>
#include <Sequence/SimData.hpp>
#include <numeric>
#include <functional>
#include <cassert>

#include <boost/container/list.hpp>
#include <boost/container/vector.hpp>
#include <boost/pool/pool_alloc.hpp>

//the type of mutation
typedef KTfwd::mutation mtype;
typedef boost::pool_allocator<mtype> mut_allocator;
typedef boost::container::list<mtype,mut_allocator > mlist;
typedef KTfwd::gamete_base<mtype,mlist> gtype;
typedef boost::pool_allocator<gtype> gam_allocator;
typedef boost::container::vector<gtype,gam_allocator > gvector;

typedef boost::unordered_set<double,boost::hash<double>,KTfwd::equal_eps > lookup_table_type;

//The mutation model returns an object of type KTfwd::mutation
KTfwd::mutation neutral_mutations_selection(gsl_rng * r,mlist * mutations,
					    const double & mu_neutral, const double & mu_selected,
					    const double & s, const double & h,
					    lookup_table_type * lookup)
{
  /*
    Choose a mutation position [0,1) that does not currently exist in the population.
  */
  double pos = gsl_rng_uniform(r);
  while( lookup->find(pos) != lookup->end() )
    {
      pos = gsl_rng_uniform(r);
    }
  lookup->insert(pos);
  assert(std::find_if(mutations->begin(),mutations->end(),boost::bind(KTfwd::mutation_at_pos(),_1,pos)) == mutations->end());
  //Choose mutation class (neutral or selected) proportional to underlying mutation rates
  bool neutral_mut = ( gsl_rng_uniform(r) <= mu_neutral/(mu_neutral+mu_selected) ) ? true : false;
  //Return a mutation of the correct type.  Neutral means s = h = 0, selected means s=s and h=h.
  return (neutral_mut == true) ? KTfwd::mutation(pos,0.,1,0.) : KTfwd::mutation(pos,s,1,h);
}

int main(int argc, char ** argv)
{
  int argument=1;
  const unsigned N = atoi(argv[argument++]);
  const double theta_neutral = atof(argv[argument++]);
  const double theta_del = atof(argv[argument++]);
  const double rho = atof(argv[argument++]);
  const double s = atof(argv[argument++]);
  const double h = atof(argv[argument++]);
  const unsigned ngens = atoi(argv[argument++]);
  const unsigned samplesize1 = atoi(argv[argument++]);
  int nreps = atoi(argv[argument++]);
  const unsigned seed = atoi(argv[argument++]);

  const double mu_neutral = theta_neutral/double(4*N);
  const double mu_del = theta_del/double(4*N);
  const double littler = rho/double(8*N);

  std::copy(argv,argv+argc,std::ostream_iterator<char*>(std::cerr," "));
  std::cerr << '\n';

  gsl_rng * r =  gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(r,seed);

  unsigned twoN = 2*N;

  while(nreps--)
    {
      //population begins monomorphic w/1 gamete and no mutations
      gvector gametes(1,gtype(twoN));
      gametes.reserve(twoN);
      mlist mutations;
      std::vector<mtype> fixations;
      std::vector<unsigned> fixation_times;
      unsigned generation;

      lookup_table_type lookup;
      for( generation = 0; generation < ngens; ++generation )
	{
	  double wbar = KTfwd::sample_diploid(r,&gametes,twoN,boost::bind(KTfwd::multiplicative_diploid(),_1,_2,2.),
					      boost::bind(KTfwd::mutation_remover(),_1,0,twoN));
	  KTfwd::remove_fixed_lost(&mutations,&fixations,&fixation_times,&lookup,generation,twoN);
	  assert(KTfwd::check_sum(gametes,twoN));

	  assert( lookup.size() == mutations.size() );
	  unsigned nmuts= KTfwd::mutate(r,&gametes,&mutations,mu_neutral+mu_del,
					boost::bind(neutral_mutations_selection,r,_1,mu_neutral,mu_del,s,h,&lookup),
					boost::bind(KTfwd::push_at_end<gtype,gvector >,_1,_2),
					boost::bind(KTfwd::insert_at_end<mtype,mlist>,_1,_2));
	  assert(KTfwd::check_sum(gametes,twoN));
	  unsigned nrec = KTfwd::recombine(r, &gametes, twoN, littler, boost::bind(gsl_rng_uniform,r));
	  assert(KTfwd::check_sum(gametes,twoN));
	}
      Sequence::SimData neutral_muts,selected_muts;

      //Take a sample of size samplesize1.  Two data blocks are returned, one for neutral mutations, and one for selected
      std::pair< std::vector< std::pair<double,std::string> >,
	std::vector< std::pair<double,std::string> > > sample = KTfwd::ms_sample_separate(r,gametes,samplesize1,twoN,true);

      neutral_muts.assign( sample.first.begin(), sample.first.end() );
      selected_muts.assign( sample.second.begin(), sample.second.end() );

      std::cout << neutral_muts << '\n' << selected_muts << '\n';
    }
}
