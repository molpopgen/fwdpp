/*
  \include bneck_selection.cc
  
  Bottleneck + exponential recovery.  With selection.  Individual-based
 */

#include <numeric>
#include <cmath>
#include <functional>
#include <cassert>
#include <iomanip>
#include <Sequence/SimData.hpp>
#include <fwdpp/diploid.hh>
//Pull mutation model from fwdpp's "sugar" layer  (@ref md_md_sugar)
#include <fwdpp/sugar/infsites.hpp>

//typedef mutation_with_age mtype;
using mtype = KTfwd::mutation;
#define SINGLEPOP_SIM
#include <common_ind.hpp>

using namespace KTfwd;

int main(int argc, char ** argv)
{
  if( argc != 14 )
    {
      std::cerr << "Error, too few arguments.\n"
		<< "Usage: " << argv[0] << ' '
		<< "N theta_neutral theta_deleterious 4Nr s h ngens N2 N3 ngens2 n nreps seed\n";
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
  const unsigned N2 = atoi(argv[argument++]);  //change N to N2 after ngens of evolution
  const unsigned N3 = atoi(argv[argument++]);  //N2 will change to N2 during ngens2 of exp. growth
  const unsigned ngens2 = atoi(argv[argument++]);
  const unsigned samplesize1 = atoi(argv[argument++]);
  int nreps = atoi(argv[argument++]);
  const unsigned seed = atoi(argv[argument++]);

  const double mu_neutral = theta_neutral/double(4*N);
  const double mu_del = theta_del/double(4*N);
  const double littler = rho/double(4*N);

  //Do some basic argument checking
  if(N2 > N)
    {
      std::cerr << "Error, N2 > N (" << N2 << " > " << N << "), but it should be N2 <= N\n";
    }
  if(N2 > N3)
    {
      std::cerr << "Error, N3 > N2 (" << N3 << " > " << N2<< "), but it should be N3 > N2\n";
    }
  if(ngens2 == 0)
    {
      std::cerr << "Error, ngens2 equals zero.  Must be > 0\n";
    }
  std::copy(argv,argv+argc,std::ostream_iterator<char*>(std::cerr," "));
  std::cerr << '\n';

  GSLrng r(seed);

  //recombination map is uniform[0,1)
  std::function<double(void)> recmap = std::bind(gsl_rng_uniform,r.get());

  unsigned twoN = 2*N;
  while(nreps--)
    {
      //Initialize a population of N diploids via KTfwd::singlepop (fwdpp/sugar/singlepop.hpp)
      singlepop_t pop(N);
      unsigned generation =  0;
      double wbar=1;
      for( generation = 0; generation < ngens; ++generation )
	{
#ifndef NDEBUG
	  for( singlepop_t::glist_t::iterator itr = pop.gametes.begin(); 
	       itr != pop.gametes.end() ; ++itr )
	    {
	      assert( itr->n > 0 );
	    }
#endif
	  assert(KTfwd::check_sum(pop.gametes,2*N));
	  wbar = KTfwd::sample_diploid(r.get(),
				       &pop.gametes,
				       &pop.diploids,
				       &pop.mutations,
				       N,
				       mu_neutral+mu_del,
				       std::bind(KTfwd::infsites(),r.get(),&pop.mut_lookup,
						 mu_neutral,mu_del,[&r](){return gsl_rng_uniform(r.get());},[&s](){return s;},[&h](){return h;}),
     				       std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,
						 std::ref(pop.neutral),std::ref(pop.selected),
						 &pop.gametes,
						 littler,
						 r.get(),
						 recmap),
				       std::bind(KTfwd::insert_at_end<singlepop_t::mutation_t,singlepop_t::mlist_t>,std::placeholders::_1,std::placeholders::_2),
				       std::bind(KTfwd::insert_at_end<singlepop_t::gamete_t,singlepop_t::glist_t>,std::placeholders::_1,std::placeholders::_2),
				       std::bind(KTfwd::multiplicative_diploid(),std::placeholders::_1,std::placeholders::_2,2.),
				       std::bind(KTfwd::mutation_remover(),std::placeholders::_1,0,2*N));
	  KTfwd::remove_fixed_lost(&pop.mutations,&pop.fixations,&pop.fixation_times,&pop.mut_lookup,generation,2*N);
	  assert(KTfwd::check_sum(pop.gametes,2*N));
	}

      //The next generation is the bottleneck
      wbar = KTfwd::sample_diploid(r.get(),
				   &pop.gametes,
				   &pop.diploids,
				   &pop.mutations,
				   N,
				   N2,
				   mu_neutral+mu_del,
				   std::bind(KTfwd::infsites(),r.get(),&pop.mut_lookup,
					     mu_neutral,mu_del,[&r](){return gsl_rng_uniform(r.get());},[&s](){return s;},[&h](){return h;}),
				   std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,
					     std::ref(pop.neutral),std::ref(pop.selected),
					     &pop.gametes,
					     littler,
					     r.get(),
					     recmap),
				   std::bind(KTfwd::insert_at_end<singlepop_t::mutation_t,singlepop_t::mlist_t>,std::placeholders::_1,std::placeholders::_2),
				   std::bind(KTfwd::insert_at_end<singlepop_t::gamete_t,singlepop_t::glist_t>,std::placeholders::_1,std::placeholders::_2),
				   std::bind(KTfwd::multiplicative_diploid(),std::placeholders::_1,std::placeholders::_2,2.),
				   std::bind(KTfwd::mutation_remover(),std::placeholders::_1,0,2*N2));
      KTfwd::remove_fixed_lost(&pop.mutations,&pop.fixations,&pop.fixation_times,&pop.mut_lookup,generation,2*N2);
      generation++;
      
      //Figure out the growth rate, etc.
      double G = std::exp( (std::log(double(N3)) - std::log(double(N2)))/double(ngens2)); 

      //Now, grow the population to size N3;
      unsigned currentN = N2,nextN;
      unsigned gens_since_bneck = 0;
      for( ; generation < ngens + ngens2 ; ++generation,++gens_since_bneck )
	{
	  nextN = round( N2*std::pow(G,gens_since_bneck+1) );
	  assert(nextN > 0);
	  wbar = KTfwd::sample_diploid(r.get(),
				       &pop.gametes,
				       &pop.diploids,
				       &pop.mutations,
				       currentN,
				       nextN,
				       mu_neutral+mu_del,
				       std::bind(KTfwd::infsites(),r.get(),&pop.mut_lookup,
						 mu_neutral,mu_del,[&r](){return gsl_rng_uniform(r.get());},[&s](){return s;},[&h](){return h;}),
     				       std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,
						 std::ref(pop.neutral),std::ref(pop.selected),
						 &pop.gametes,
						 littler,
						 r.get(),
						 recmap),
				       std::bind(KTfwd::insert_at_end<singlepop_t::mutation_t,singlepop_t::mlist_t>,std::placeholders::_1,std::placeholders::_2),
				       std::bind(KTfwd::insert_at_end<singlepop_t::gamete_t,singlepop_t::glist_t>,std::placeholders::_1,std::placeholders::_2),
				       std::bind(KTfwd::multiplicative_diploid(),std::placeholders::_1,std::placeholders::_2,2.),
				       std::bind(KTfwd::mutation_remover(),std::placeholders::_1,0,2*nextN));
	  KTfwd::remove_fixed_lost(&pop.mutations,&pop.fixations,&pop.fixation_times,&pop.mut_lookup,generation,2*nextN);
	  currentN=nextN;
	}
      Sequence::SimData neutral_muts,selected_muts;
      
      //Take a sample of size samplesize1.  Two data blocks are returned, one for neutral mutations, and one for selected
      std::pair< std::vector< std::pair<double,std::string> >,
		 std::vector< std::pair<double,std::string> > > sample = ms_sample_separate(r.get(),&pop.diploids,samplesize1);
      
      neutral_muts.assign( sample.first.begin(), sample.first.end() );
      selected_muts.assign( sample.second.begin(), sample.second.end() );
      
      std::cout << neutral_muts << '\n' << selected_muts << '\n';
    }
}
