/*! \include RHH.cc
  Recurrent hitch-hiking model of Braverman et al, extended to allow selected mutations within the region
 */

#include <fwdpp/diploid.hh>
#include <Sequence/SimData.hpp>
#include <boost/unordered_set.hpp>

#include <boost/container/list.hpp>
#include <boost/container/vector.hpp>
#include <boost/pool/pool_alloc.hpp>

struct mutation_with_origin : public KTfwd::mutation
//records the generation when the mutation entered the population
{
  mutable unsigned o;
  mutation_with_origin (const unsigned & __o ,const double & position, const double & sel_coeff,const unsigned & count,
			const double & dominance = 0.5) :
    o(__o),KTfwd::mutation(position,sel_coeff,count,dominance)
  {
  }
};

typedef mutation_with_origin mtype;
typedef boost::pool_allocator<mtype> mut_allocator;
typedef boost::container::list<mtype,mut_allocator > mlist;
typedef KTfwd::gamete_base<mtype,mlist> gtype;
typedef boost::pool_allocator<gtype> gam_allocator;
typedef boost::container::vector<gtype,gam_allocator > gvector;
typedef boost::unordered_set<double,boost::hash<double>,KTfwd::equal_eps > lookup_table_type;


mtype RHH_mutation_model( gsl_rng * r, const unsigned & generation, const double mu_neutral, const double & mu_selected,
			  const double & p_selected_within,
			  const double & s, const double & h, const double & maxd, lookup_table_type * lookup ) //keep h (dominance) to 1 to match theory
{
  //Will the mutation be neutral, or a selected one?
  double neutral = (gsl_rng_uniform(r) <= mu_neutral/(mu_neutral+mu_selected)) ? true : false;
  double pos = std::numeric_limits<double>::max(); //set to some impossible dummy value
  double smut=0.;
  if( neutral )
    {
      pos = gsl_rng_uniform(r);
      while( lookup->find(pos) != lookup->end() )
	{
	  pos = gsl_rng_uniform(r);
	}
    }
  else
    {
      bool within = (gsl_rng_uniform(r) <= p_selected_within) ? true : false;
      pos = (within) ? gsl_rng_uniform(r) : gsl_ran_flat(r,0.,maxd);
      if( ! within )
	{
	  if (gsl_rng_uniform(r) <= 0.5 )
	    {
	      pos *= -1.;
	    }
	  else
	    {
	      pos += 1.;
	    }
	}
      while( lookup->find(pos) != lookup->end() )
	{
	  pos = (within) ? gsl_rng_uniform(r) : gsl_ran_flat(r,0.,maxd);
	  if( ! within )
	    {
	      if (gsl_rng_uniform(r) <= 0.5)
		{
		  pos *= -1.;
		}
	      else
		{
		  pos += 1.;
		}
	    }
	}
      smut=s;
    }
  lookup->insert(pos);
  return mtype(generation,pos,smut,1,h);
}

double recurrent_sweep_genetic_map(gsl_rng * r, const double & littler_neut,
				   const double & ttl_rec_rate,
				   const double & maxd)
{
  double rdm = gsl_rng_uniform(r),pos;
  if( rdm <= littler_neut/ttl_rec_rate )
    {
      pos = gsl_rng_uniform(r);
    }
  else
    {
      pos = (gsl_rng_uniform(r) <= 0.5) ? gsl_ran_flat(r,-1.*maxd,0.) : gsl_ran_flat(r,1.,1.+maxd);
    }
  return pos;
}
 
int main(int argc, char ** argv)
{
  if ( argc != 12 )
    {
      std::cerr << "Too few arguments.\n"
		<< "Usage: RHH N theta rho nsites s Lambda burnin ngens samplesize nreps seed\n";
      exit(10);
    }
  int argument=1;
  const unsigned N = atoi(argv[argument++]);           //Number of diploids
  const double theta = atof(argv[argument++]);         //4*n*mutation rate for neutral region.  Note: mutation rate is per REGION, not SITE!!
  const double rho = atof(argv[argument++]);           //4*n*recombination rate for neutral region.  Note: recombination rate is per REGION, not SITE!!
  const unsigned nsites = atoi(argv[argument++]);               //length of neutral region in bp
  const double s = atof(argv[argument++]);                 //selection coefficient
  const double h = 1;                                     //codominance w/fitnesses 1,1+s,1+2s
  const double lambda = atof(argv[argument++]);            //4NLambda = E[# sweeps] per 4N gens per site
  const unsigned burnin = atoi(argv[argument++]);       //Number of generations to simulate w/o selection
  const unsigned ngens = atoi(argv[argument++]);       //Number of generations to simulate w/selection
  const unsigned samplesize1 = atoi(argv[argument++]); //Sample size to draw from the population
  int nreps = atoi(argv[argument++]);                  //Number of replicates to simulate
  const unsigned seed = atoi(argv[argument++]);        //Random number seed

  //4Ns is the max (genetic) distance where a sweep may have an effect on a neutrla site (Kaplan et al, and others)
  const double maxd = 4.*double(N)*s;
  //neutral mutation rate (per gamete per gen)
  const double mu_neutral = theta/double(4*N);
  //work backwards from formula for fixation prob. under genic selection to mutation rate to beneficial alleles
  const double mu_pos_bp = (lambda > 0.) ? ((lambda*(1-exp(-4.*double(N)*s)))/(8.*pow(double(N),2.)*(1.-exp(-2.*s)))) : 0.;

  //mutation rate to selected sites w/in region where neutral mutations occur
  const double mu_pos_in = mu_pos_bp*double(nsites);
  //maxdbp maxd in base pairs, using length of sampled region to scale genetic to physical distance
  const double maxdbp = maxd/(rho/double(nsites-1));
  //total mutation rate to selected muts
  const double mu_selected_ttl =mu_pos_in + mu_pos_bp*2.*maxdbp; 
  
  //recombination rates
  const double littler_neut = (rho/double(8*N)); //within sampled region
  //on entire part of chromo being considered
  const double littler = littler_neut + 2.*littler_neut*(double(maxdbp)/double(nsites-1));

  std::cerr << '#';
  std::copy(argv,argv+argc,std::ostream_iterator<char*>(std::cerr," "));
  std::cerr << '\n';

  //Initiate random number generation system
  gsl_rng * r =  gsl_rng_alloc(gsl_rng_ranlxs2);
  gsl_rng_set(r,seed);

  unsigned twoN = 2*N;

  std::cerr.precision(10.);

  while(nreps--)
    {
      //the population begins with 1 gamete with no mutations amd initial count 2N
      gvector gametes(1,gtype(twoN));
      mlist mutations;  //the population is devoid of mutations
      std::vector<mtype> fixations;  //store mutation that fix.  Passed to KTfwd::remove_fixed_lost
      std::vector<unsigned> fixation_times; //store when they fix.  Passed to KTfwd::remove_fixed_lost
      unsigned generation;
      double wbar;
      lookup_table_type lookup;  //this is our lookup table for the mutation model

      for( generation = 0; generation < burnin; ++generation )
	{
	  wbar = KTfwd::sample_diploid(r,&gametes,twoN,
				       std::bind(KTfwd::multiplicative_diploid(),std::placeholders::_1,std::placeholders::_2,2.),
				       std::bind(KTfwd::mutation_remover(),std::placeholders::_1,0,twoN));       
	  KTfwd::remove_fixed_lost(&mutations,&fixations,&fixation_times,&lookup,generation,twoN);
	  assert( lookup.size() == mutations.size() );
	  assert(KTfwd::check_sum(gametes,twoN));
	  unsigned nmuts = KTfwd::mutate(r,&gametes,&mutations,mu_neutral,
					 std::bind(RHH_mutation_model,r,generation,mu_neutral,0.,0.,0.,h,maxd,&lookup),
					 std::bind(KTfwd::push_at_end<gtype,gvector >,std::placeholders::_1,std::placeholders::_2),
					 std::bind(KTfwd::insert_at_end<mtype,mlist>,std::placeholders::_1,std::placeholders::_2));
	  assert( lookup.size() == mutations.size() );
	  assert(KTfwd::check_sum(gametes,twoN));
	  unsigned nrec = KTfwd::recombine(r, 
					   &gametes,
					   twoN, 
					   littler_neut, 
					   std::bind(gsl_rng_uniform,r));
	  assert(KTfwd::check_sum(gametes,twoN));
	}
      for(  ;generation < ngens+burnin; ++generation )
	{
	  wbar = KTfwd::sample_diploid(r,&gametes,twoN,
				       std::bind(KTfwd::multiplicative_diploid(),std::placeholders::_1,std::placeholders::_2,2.),
				       std::bind(KTfwd::mutation_remover(),std::placeholders::_1,0,twoN));       
	  KTfwd::remove_fixed_lost(&mutations,&fixations,&fixation_times,&lookup,generation,twoN);
	  assert( lookup.size() == mutations.size() );
	  assert(KTfwd::check_sum(gametes,twoN));
	  unsigned nmuts = KTfwd::mutate(r,&gametes,&mutations,mu_neutral+mu_selected_ttl,
					 std::bind(RHH_mutation_model,r,generation,mu_neutral,mu_selected_ttl,
						     mu_pos_in/mu_selected_ttl,s,h,maxd,&lookup),
					 std::bind(KTfwd::push_at_end<gtype,gvector >,std::placeholders::_1,std::placeholders::_2),
					 std::bind(KTfwd::insert_at_end<mtype,mlist>,std::placeholders::_1,std::placeholders::_2));
	  assert( lookup.size() == mutations.size() );
	  assert(KTfwd::check_sum(gametes,twoN));
	  //only recombine in neutral region unless selected mutations are present...
	  bool selected=false;
	  for(mlist::iterator itr = mutations.begin();selected == false && itr!=mutations.end() ; ++itr )
	    {
	      selected = !itr->neutral;
	    }
	  if(selected)
	    {
	      unsigned nrec = KTfwd::recombine(r, 
					       &gametes,
					       twoN, 
					       littler,
					       //genetic map uniform over total region, neutral + selected
					       std::bind( recurrent_sweep_genetic_map,r,littler_neut,
							    littler,maxd ) );
	    }
	  else
	    {
	      unsigned nrec = KTfwd::recombine(r, 
					       &gametes,
					       twoN, 
					       littler_neut, 
					       //genetic map uniform over sampled region
					       std::bind(gsl_rng_uniform,r));
	    }
	  assert(KTfwd::check_sum(gametes,twoN));
	}
      Sequence::SimData sdata;
      //.first is the neutral data, .second the selected sites
      std::pair< std::vector<std::pair<double,std::string> >,
	std::vector<std::pair<double,std::string> > > mslike = KTfwd::ms_sample_separate(r,gametes,samplesize1,twoN,true);

      //Convert the results of KTfwd::ms_sample to a Sequence::SimData object and print to the screen
      if(!mslike.first.empty())
	{
	  //only output the neutral sites, to mimic coalescent sims
	  //where we look @ neutral variability only to study the effect of linked selection
	  sdata.assign(mslike.first.begin(),mslike.first.end());
	  std::cout << sdata << '\n';
	}
      else
	{
	  std::cout << "//\nsegsites: 0\n";
	}  
      //write the origination and fixation times of all non-neutral substitutions to STDERR
      for(unsigned i = 0 ; i < fixations.size() ; ++i)
	{
	  if( ! fixations[i].neutral )
	    {
	      std::cerr << fixations[i].pos << '\t' << fixations[i].o << '\t' << fixation_times[i] << '\n';
	    }
	}
    }
}
