/*
  Pseudo separate sex model with sex-specific mutation rates 
  and fitness effects of mutations.

  NOTE: this is a toy example.  It is possible for two "males"
  or two "females" to be chosen as parents.  The point here is 
  to show how custom diploid types may be used, and how that extra
  data may be used to do things like calculate fitness and (diploid label) 
  by (mutation label) interactions, etc.

  \include sex_limited.cc
  Total madness:
  1. Opposite fitness effects in "males" vs. "females"
  2. Gaussian stabiliziing selection
  3. House of cards
  Interaction parameters are hard-coded in.
*/

#include <limits>
#include <algorithm>

#include <boost/container/list.hpp>
#include <boost/container/vector.hpp>
#include <boost/pool/pool_alloc.hpp>
#include <boost/unordered_set.hpp>
#include <boost/functional/hash.hpp>

//Main fwdpp library header
#include <fwdpp/diploid.hh>
//Include the necessary "sugar" components
//We need to get the 'deep' version of singlepop, as we need to make a custom singlepop_serialized_t for our sim
#include <fwdpp/sugar/singlepop/singlepop.hpp>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/sugar/serialization.hpp>
#include <fwdpp/sugar/infsites.hpp>
#include <fwdpp/experimental/sample_diploid.hpp>
//FWDPP-related stuff

struct sex_specific_mutation : public KTfwd::mutation_base
{
  //Effect size
  double s;
  //The effect size applies to this sex, is 0 otherwise
  bool sex;
  sex_specific_mutation(const double & __pos, const double & __s, const bool & __sex,
			const unsigned & __n, const bool & __neutral)
    : KTfwd::mutation_base(__pos,__n,__neutral),s(__s),sex(__sex)
  {	
  }
};

using mtype = sex_specific_mutation;
using mlist_t = boost::container::list<mtype,boost::pool_allocator<mtype> >;
using gamete_t = KTfwd::gamete_base<mtype,mlist_t>;
using glist_t = boost::container::list<gamete_t, boost::pool_allocator<gamete_t>>;

//We need to define a custom diploid genotype for our model
struct diploid_t : public KTfwd::tags::custom_diploid_t
{
  using first_type = glist_t::iterator;
  using second_type = glist_t::iterator;
  first_type first;
  second_type second;
  bool sex; //male = false, I guess...
  //constructors, etc.
  diploid_t() : first(first_type()),second(second_type()),sex(false) {}
  //"perfect forwarding" constructor does not work with iterator from boost containers...
  //diploid_t(first_type && g1, first_type && g2) : first(std::forward(g1)),second(std::forward(g2)),i(numeric_limits<unsigned>::max()) {}
  diploid_t(first_type g1, first_type g2) : first(g1),second(g2),sex(false) {}
  //The following constructors SHOULD be generated automagically by your compiler, so you don't have to:
  //(no idea what, if any, performance effect this may have.  Worst case is prob. the move constructor doesn't get auto-generated...
  //diploid_t( const diploid_t & ) = default;
  //diploid_t( diploid_t && ) = default;
  //diploid_t & operator=(const diploid_t &) = default;
};

//Define our our population type via KTfwd::sugar 
using poptype = KTfwd::sugar::singlepop_serialized<mtype,KTfwd::mutation_writer,KTfwd::mutation_reader<mtype>,
						   mlist_t,
						   glist_t,
						   boost::container::vector< diploid_t >,
						   boost::container::vector<mtype>,
						   boost::container::vector<unsigned>,
						   boost::unordered_set<double,boost::hash<double>,KTfwd::equal_eps>
						   >;

/*
  We will use a gsl_rng_mt19937 as our RNG.
  This type is implicitly convertible to gsl_rng *,
  and auto-handles the gsl_rng_free steps, etc.
*/
using GSLrng = KTfwd::GSLrng_t<KTfwd::GSL_RNG_MT19937>;

mtype sex_specific_mut_model( gsl_rng * r,
			      poptype::mlist_t * mutations,
			      poptype::lookup_table_t * lookup,
			      const double & mu_total,
			      const double & mu_male,
			      const double & mu_female,
			      const double & sigma )
{
  double pos = gsl_rng_uniform(r);
  while(lookup->find(pos) != lookup->end())
    {
      pos = gsl_rng_uniform(r);
    }
  lookup->insert(pos);
  double u = gsl_rng_uniform(r);
  if(u <= mu_male/mu_total)
    {
      return mtype(pos,gsl_ran_gaussian(r,sigma),false,1,false);
    }
  else if (u <= (mu_male+mu_female)/mu_total)
    {
      return mtype(pos,gsl_ran_gaussian(r,sigma),true,1,false);
    }
  //Otherwise, neutral mutation
  //We "hack" this and assign the mutation a "male" type,
  //As they'll never be used in a fitness calc,
  //as they'll be stored in mutations rather than
  //smutations
  return mtype(pos,0.,false,1,true);
}

//We need a fitness model
double sex_specific_fitness( const poptype::dipvector_t::const_iterator & dip, gsl_rng * r, const double & sigmaE )
{
  double trait_value = std::accumulate( dip->first->smutations.begin(),
					dip->first->smutations.end(),
					0.,[&dip](const double & a, const poptype::mlist_t::const_iterator & m)
					{
					  return a + ((dip->sex==m->sex) ? m->s : 0.);
					} );
  trait_value += std::accumulate( dip->second->smutations.begin(),
				  dip->second->smutations.end(),
				  0.,[&dip](const double & a, const poptype::mlist_t::const_iterator & m)
				  {
				    return a + ((dip->sex==m->sex) ? m->s : 0.);
				  } );
  return std::exp( -std::pow(trait_value+gsl_ran_gaussian(r,sigmaE),2.)/2.);
}


/*! \brief Continuous Snowdrift Game from Doebeli, Hauert, and Killingback (2004, Science, 306:859--862)
  as implemented by Wakano and Lehmann (2014, J Theor Biol, 351:83--95)
  \param phenotypes phenotypes in population
  \param fitnesses fitnesses in population
  \param b1 linear benefit term
  \param b2 quadratic benefit term
  \param c1 linear cost term
  \param c2 quadratic cost term
  \return Mean payoff (fitness) from pairwise interactions in snowdrift game
  \ingroup fitness
*/

int main(int argc, char ** argv)
{
  if (argc != 12)
    {
      std::cerr << "Too few arguments.\n"
		<< "Usage: " << argv[0]
		<< " N mu_neutral mu_male mu_female sigma_mu sigma_e recrate ngens samplesize nreps seed\n";
      exit(10);
    }
  int argument=1;
  
  const unsigned N = atoi(argv[argument++]);
  const double mu_neutral = atof(argv[argument++]);
  const double mu_male = atof(argv[argument++]);
  const double mu_female = atof(argv[argument++]);
  const double sigma = atof(argv[argument++]);
  const double sigmaE = atof(argv[argument++]);
  const double recrate = atof(argv[argument++]);
  const unsigned ngens = atoi(argv[argument++]);
  const unsigned samplesize1 = atoi(argv[argument++]);
  const unsigned nreps = atoi(argv[argument++]);
  const unsigned seed = atoi(argv[argument++]);

  GSLrng rng(seed);
  std::function<double(void)> recmap = std::bind(gsl_rng_uniform,rng); //uniform crossover map
  const double mu_total = mu_neutral+mu_male+mu_female;
  for( unsigned rep = 0 ; rep < nreps ; ++rep )
    {
      poptype pop(N);
      for( unsigned generation = 0 ; generation < ngens ; ++generation )
	{
	  //Assign "sex"
	  for( auto dip = pop.diploids.begin() ; dip != pop.diploids.end() ; ++dip )
	    {
	      dip->sex = (gsl_rng_uniform(rng) <= 0.5); //false = male, true = female.
	    }
	  //double wbar = KTfwd::sample_diploid(rng,
	  double wbar = KTfwd::experimental::sample_diploid(rng,
							    &pop.gametes,
							    &pop.diploids,
							    &pop.mutations,
							    N,
							    mu_total,
							    std::bind(sex_specific_mut_model,rng.r.get(),std::placeholders::_1,
								      &pop.mut_lookup,mu_total,mu_male,mu_female,sigma),
							    std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,
								      &pop.gametes,
								      recrate, 
								      rng,
								      recmap),
							    std::bind(KTfwd::insert_at_end<poptype::mutation_t,poptype::mlist_t>,std::placeholders::_1,std::placeholders::_2),
							    std::bind(KTfwd::insert_at_end<poptype::gamete_t,poptype::glist_t>,std::placeholders::_1,std::placeholders::_2),
							    std::bind(sex_specific_fitness,std::placeholders::_1,rng,sigmaE),
							    std::bind(KTfwd::mutation_remover(),std::placeholders::_1,0,2*pop.N));
	  KTfwd::remove_fixed_lost(&pop.mutations,&pop.fixations,&pop.fixation_times,&pop.mut_lookup,generation,2*pop.N);
	}
      Sequence::SimData neutral_muts,selected_muts;
      
      //Take a sample of size samplesize1.  Two data blocks are returned, one for neutral mutations, and one for selected
      std::pair< std::vector< std::pair<double,std::string> >,
		 std::vector< std::pair<double,std::string> > > sample = ms_sample_separate(rng,&pop.diploids,samplesize1);
      
      neutral_muts.assign( sample.first.begin(), sample.first.end() );
      selected_muts.assign( sample.second.begin(), sample.second.end() );

      std::cout << neutral_muts << '\n' << selected_muts << '\n';
    }
}
