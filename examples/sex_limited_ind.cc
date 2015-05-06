/*
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

//FWDPP-related stuff
using mtype = KTfwd::mutation;
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
  if (argc != 11)
    {
      std::cerr << "Too few arguments.\n"
		<< "Usage: " << argv[0]
		<< " N mu_neutral mu_deleterious recrate s h ngens samplesize nreps seed\n";
      exit(10);
    }
  int argument=1;
  
  const unsigned N = atoi(argv[argument++]);
  const double mu_neutral = atof(argv[argument++]);
  const double mu_del = atof(argv[argument++]);
  const double recrate = atof(argv[argument++]);
  const double s = atof(argv[argument++]);
  const double h = atof(argv[argument++]);
  const unsigned ngens = atoi(argv[argument++]);
  const unsigned samplesize1 = atoi(argv[argument++]);
  const unsigned nreps = atoi(argv[argument++]);
  const unsigned seed = atoi(argv[argument++]);

  GSLrng rng(seed);
  std::function<double(void)> recmap = std::bind(gsl_rng_uniform,rng); //uniform crossover map
  for( unsigned rep = 0 ; rep < nreps ; ++rep )
    {
      poptype pop(N);
      for( unsigned generation = 0 ; generation < ngens ; ++generation )
	{
	  //Fill phenotypes
	  unsigned i = 0;
	  double wbar = KTfwd::sample_diploid(rng,
					      &pop.gametes,
					      &pop.diploids,
					      &pop.mutations,
					      N,
					      mu_neutral+mu_del,
					      std::bind(KTfwd::infsites(),rng,std::placeholders::_1,&pop.mut_lookup,
							mu_neutral,mu_del,[&rng](){return gsl_rng_uniform(rng);},
							[&rng,&s](){return -1.*gsl_ran_exponential(rng,s);},
							[&h](){return h;}),
					      std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,
							&pop.gametes,
							recrate, 
							rng,
							recmap),
					      std::bind(KTfwd::insert_at_end<poptype::mutation_t,poptype::mlist_t>,std::placeholders::_1,std::placeholders::_2),
					      std::bind(KTfwd::insert_at_end<poptype::gamete_t,poptype::glist_t>,std::placeholders::_1,std::placeholders::_2),
					      //std::bind(snowdrift_diploid(),std::placeholders::_1,std::cref(phenotypes),b1,b2,c1,c2),
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
