/*
  \include diploid_fixed_sh_ind.cc
  
  Same as diploid_fixed_sh.cc, but individual-based
 */

#include <fwdpp/diploid.hh>
#include <Sequence/SimData.hpp>
#include <numeric>
#include <functional>
#include <cassert>
#include <iomanip>

//the type of mutation
typedef KTfwd::mutation mtype;
#include <common_ind.hpp>

using namespace KTfwd;

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
  assert(std::find_if(mutations->begin(),mutations->end(),std::bind(KTfwd::mutation_at_pos(),std::placeholders::_1,pos)) == mutations->end());
  //Choose mutation class (neutral or selected) proportional to underlying mutation rates
  bool neutral_mut = ( gsl_rng_uniform(r) <= mu_neutral/(mu_neutral+mu_selected) ) ? true : false;
  //Return a mutation of the correct type.  Neutral means s = h = 0, selected means s=s and h=h.
  return (neutral_mut == true) ? KTfwd::mutation(pos,0.,1,0.) : KTfwd::mutation(pos,s,1,h);
}

int main(int argc, char ** argv)
{
  if (argc != 11)
    {
      std::cerr << "Too few arguments.\n"
		<< "Usage: diploid_fixed_sh_ind_lambda N theta_neutral theta_deleterious rho s h ngens samplesize nreps seed\n";
      exit(10);
    }
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
  const double littler = rho/double(4*N);

  std::copy(argv,argv+argc,std::ostream_iterator<char*>(std::cerr," "));
  std::cerr << '\n';

  gsl_rng * r =  gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(r,seed);

  unsigned twoN = 2*N;

  KTfwd::site_dependent_fitness fitness_model;
  while(nreps--)
    {
      //population begins monomorphic w/1 gamete and no mutations
      glist gametes(1,gtype(twoN));
      std::vector< std::pair< glist::iterator, glist::iterator > > diploids(N,
									    std::make_pair(gametes.begin(),
											   gametes.begin()));
      mlist mutations;
      std::vector<mtype> fixations;
      std::vector<unsigned> fixation_times;
      unsigned generation;

      lookup_table_type lookup;
      double wbar=1;
      for( generation = 0; generation < ngens; ++generation )
	{
#ifndef NDEBUG
	  for( glist::iterator itr = gametes.begin(); 
	       itr != gametes.end() ; ++itr )
	    {
	      assert( itr->n > 0 );
	    }
#endif
	  assert(KTfwd::check_sum(gametes,2*N));
	  wbar = KTfwd::sample_diploid(r,
				       &gametes,
				       &diploids,
				       &mutations,
				       N,
				       mu_neutral+mu_del,
				       //The mutation model function will be passed a non-const pointer to an mlist
				       [&](mlist * __mutations){ return neutral_mutations_selection(r,__mutations,mu_neutral,mu_del,s,h,&lookup);},
				       //The recombination policy must take two non-const iterators from the glist
				       [&](glist::iterator & g1,
					   glist::iterator & g2) { return KTfwd::recombine_gametes(r,littler,&gametes,g1,g2,
												   //This nested lambda is our genetic map: uniform on interval (0,1]
												   [&](){return gsl_rng_uniform(r);}); },
					   //glist::iterator & g2) { return KTfwd::recombine_gametes(r,littler,&gametes,g1,g2,recmap); },
				       //The mutation insertion policy takes a const mtype and a non-const pointer to an mlist
				       [](const mtype & m,mlist * __mutations) { return __mutations->insert(__mutations->end(),m); },
				       //The gamete insertion policy takes a const gtype and a non-const pointer to a glist
				       []( const gtype & g,glist * __gametes) { return __gametes->insert(__gametes->end(),g); },
				       //Our fitness model takes two const glist::const_iterators  and two policies as arguments
				       [&](const glist::const_iterator & g1, const glist::const_iterator &g2) {
					 return std::cref(fitness_model)(g1,g2,
									 //The first policy is what to do for a homozygous mutation
									 //The "homozygote" policy must take a non-const reference to a double
									 //and a const mlist::iterator as arguments
									 [](double & fitness,const mlist::iterator & mut) { 
									   fitness *= (1. + 2.*mut->s);
									 },
									 //The second is what to do for a heterozygous mutation
									 //The "heterozygote" policy must take a non-const reference to a double
									 //and a const mlist::iterator as arguments
									 [](double & fitness,const mlist::iterator & mut) {
									   fitness *= (1. + mut->h*mut->s);
									 });},
				       //Finally, mutations will be removed from gametes if they are extinct or fixed
				       [&N](const mlist::iterator & i){ return (!i->n || i->n == 2*N); });
	  KTfwd::remove_fixed_lost(&mutations,&fixations,&fixation_times,&lookup,generation,2*N);
	  assert(KTfwd::check_sum(gametes,2*N));
	}
      Sequence::SimData neutral_muts,selected_muts;

      //Take a sample of size samplesize1.  Two data blocks are returned, one for neutral mutations, and one for selected
      std::pair< std::vector< std::pair<double,std::string> >,
	std::vector< std::pair<double,std::string> > > sample = ms_sample_separate(r,&diploids,samplesize1);

      neutral_muts.assign( sample.first.begin(), sample.first.end() );
      selected_muts.assign( sample.second.begin(), sample.second.end() );

      std::cout << neutral_muts << '\n' << selected_muts << '\n';
    }
}
