/*! \include diploid_twopop_mig.cc
  Population split followed by migration.
 */
#include <config.h>
#include <fwdpp/diploid.hh>

#include <numeric>
#include <functional>
#include <cassert>

#include <Sequence/SimData.hpp>

struct mutation_with_age : public KTfwd::mutation_base
{
  unsigned g;
  double s,h;
  mutation_with_age(const unsigned & __generation,const double & position, const unsigned & count, const bool & isneutral = true)
    : KTfwd::mutation_base(position,count,isneutral),g(__generation),s(0.),h(0.)
  {	
  }
};

typedef mutation_with_age mtype;
#include <common_gamete.hpp>

#if defined(HAVE_BOOST_VECTOR) && defined(HAVE_BOOST_LIST) && defined(HAVE_BOOST_UNORDERED_SET) && defined(HAVE_BOOST_POOL_ALLOC) && defined(HAVE_BOOST_HASH) && !defined(USE_STANDARD_CONTAINERS)
typedef boost::container::vector< gvector > mpop_container;
#else
typedef std::vector< gvector > mpop_container;
#endif

mutation_with_age neutral_mutations_inf_sites(gsl_rng * r,const unsigned & generation,mlist * mutations,
					      lookup_table_type * lookup)
{
  double pos = gsl_rng_uniform(r);
  while(lookup->find(pos) != lookup->end())
    {
      pos = gsl_rng_uniform(r);
    }
  lookup->insert(pos);
  assert(std::find_if(mutations->begin(),mutations->end(),std::bind(KTfwd::mutation_at_pos(),std::placeholders::_1,pos)) == mutations->end());
  return mutation_with_age(generation,pos,1,true);
}

Sequence::SimData merge( const std::vector<std::pair<double,std::string> > & sample1,
			 const std::vector<std::pair<double,std::string> > & sample2,
			 const unsigned & nsam );

int main(int argc, char ** argv)
{
  if (argc != 10)
    {
      std::cerr << "Too few arguments.\n"
		<< "Usage: diploid_twopop_mig N theta rho ngens_b4_split ngens_after_split M samplesize nreps seed\n";
      exit(10);
    }
  int argument=1;
  const unsigned N = atoi(argv[argument++]);
  const double theta = atof(argv[argument++]);
  const double rho = atof(argv[argument++]);
  const unsigned ngens1 = atoi(argv[argument++]); //generations to evolve before split
  const unsigned ngens2 = atoi(argv[argument++]); //generations to evolve after split
  const double M = atof(argv[argument++]); //4Nm -- symmetric b/w pops
  const unsigned samplesize1 = atoi(argv[argument++]);
  int nreps = atoi(argv[argument++]);
  const unsigned seed = atoi(argv[argument++]);

  const double mu = theta/double(4*N);
  const double littler = rho/double(8*N);
  const double littlem = M/double(4*N);

  gsl_rng * r =  gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(r,seed);

  unsigned twoN = 2*N;

  //create a vector of fitness functions for each population
  std::vector<std::function<double (gvector::const_iterator,
				      gvector::const_iterator)> > vbf;
  vbf.push_back(std::bind(KTfwd::multiplicative_diploid(),std::placeholders::_1,std::placeholders::_2,2.));
  vbf.push_back(std::bind(KTfwd::multiplicative_diploid(),std::placeholders::_1,std::placeholders::_2,2.));

  while(nreps--)
    {
      //the population begins with 1 gamete with no mutations
      mpop_container metapop(1,gvector(1,gtype(twoN)));
      mlist mutations;
      std::vector<mtype> fixations;
      std::vector<unsigned> fixation_times;
      unsigned generation;

      double wbar;
      lookup_table_type lookup;
      for( generation = 0; generation < ngens1; ++generation )
	{

	  wbar = KTfwd::sample_diploid(r,&metapop[0],twoN,std::bind(KTfwd::multiplicative_diploid(),std::placeholders::_1,std::placeholders::_2,2.),
				   std::bind(KTfwd::mutation_remover(),std::placeholders::_1,0,twoN));
	  KTfwd::remove_fixed_lost(&mutations,&fixations,&fixation_times,&lookup,generation,twoN);
	  assert(KTfwd::check_sum(metapop[0],twoN));

	  unsigned nmuts= KTfwd::mutate(r,&metapop[0],&mutations,mu,
	  				std::bind(neutral_mutations_inf_sites,r,generation,std::placeholders::_1,&lookup),
					std::bind(KTfwd::push_at_end<gtype,gvector >,std::placeholders::_1,std::placeholders::_2),
					std::bind(KTfwd::insert_at_end<mtype,mlist>,std::placeholders::_1,std::placeholders::_2));

	  assert(KTfwd::check_sum(metapop[0],twoN));
	  unsigned nrec = KTfwd::recombine(r, &metapop[0], twoN, littler, std::bind(gsl_rng_uniform,r));
	  assert(KTfwd::check_sum(metapop[0],twoN));
	}

      //make another copy of population 0 and add it to metapop
      //NOTE: with boost containers, for some reason, we need explicit constructor call in the push_back:
      metapop.push_back(gvector(metapop[0]));
      std::vector<unsigned> twoNs(2,twoN);
      unsigned ttlsize = std::accumulate(twoNs.begin(),twoNs.end(),0u);
      assert(ttlsize == 2*twoN);
      //start sampling the metapopulation
      for( ; generation < ngens1 + ngens2 ; ++generation )
	{
	  //migrate
	  KTfwd::migrate(r,&metapop[0],&metapop[1],twoN,twoN,littlem);
	  std::vector<double> wbars = KTfwd::sample_diploid(r,&metapop,&twoNs[0],ttlsize,
							    vbf,   //vector of fitness functions
							    std::bind(KTfwd::mutation_remover(),std::placeholders::_1,0,twoN));
	  assert(KTfwd::check_sum(metapop[0],twoN));
	  assert(KTfwd::check_sum(metapop[1],twoN));
	  KTfwd::remove_fixed_lost(&mutations,&fixations,&fixation_times,generation,2*twoN);

	  unsigned nmuts= KTfwd::mutate(r,&metapop[0],&mutations,mu,
	  				std::bind(neutral_mutations_inf_sites,r,generation,std::placeholders::_1,&lookup),
					std::bind(KTfwd::push_at_end<gtype,gvector >,std::placeholders::_1,std::placeholders::_2),
					std::bind(KTfwd::insert_at_end<mtype,mlist>,std::placeholders::_1,std::placeholders::_2));
	  
	  assert(KTfwd::check_sum(metapop[0],twoN));
	  unsigned nrec = KTfwd::recombine(r, &metapop[0], twoN, littler, std::bind(gsl_rng_uniform,r));
	  assert(KTfwd::check_sum(metapop[0],twoN));
	  
	  nmuts= KTfwd::mutate(r,&metapop[1],&mutations,mu,
			       std::bind(neutral_mutations_inf_sites,r,generation,std::placeholders::_1,&lookup),
			       std::bind(KTfwd::push_at_end<gtype,gvector >,std::placeholders::_1,std::placeholders::_2),
			       std::bind(KTfwd::insert_at_end<mtype,mlist>,std::placeholders::_1,std::placeholders::_2));
	  
	  assert(KTfwd::check_sum(metapop[1],twoN));
	  nrec = KTfwd::recombine(r, &metapop[1], twoN, littler, std::bind(gsl_rng_uniform,r));
	  assert(KTfwd::check_sum(metapop[1],twoN));

	} 
      
      //sample each pop, then merge the results into 1 "ms" block:
      std::vector<std::pair<double,std::string> > sample1 = KTfwd::ms_sample(r,metapop[0],samplesize1,twoN,false);
      std::vector<std::pair<double,std::string> > sample2 = KTfwd::ms_sample(r,metapop[1],samplesize1,twoN,false);
      Sequence::SimData tsample=merge(sample1,sample2,samplesize1);
      std::cout << tsample << '\n';
    }
}

Sequence::SimData merge( const std::vector<std::pair<double,std::string> > & sample1,
			 const std::vector<std::pair<double,std::string> > & sample2 ,
			 const unsigned & nsam)
{
  std::map<double, std::string> temp;

  for( unsigned i=0;i<sample1.size();++i)
    {
      temp[sample1[i].first] = std::string(sample1[i].second + std::string(nsam,'0'));
    }

  for( unsigned i=0;i<sample2.size();++i)
    {
      std::map<double,std::string>::iterator itr = temp.find(sample2[i].first);
      if( itr == temp.end() )
	{
	  temp[sample2[i].first] = std::string( std::string(nsam,'0') + sample2[i].second );
	}
      else
	{
	  std::copy( sample2[i].second.begin(),sample2[i].second.end(),itr->second.begin()+nsam );
	}
    }
  std::vector<std::pair<double,std::string> > rv( temp.begin(),temp.end() );
  std::sort(rv.begin(),rv.end(),
		[](std::pair<double,std::string> lhs,
		   std::pair<double,std::string> rhs) { return lhs.first < rhs.first; });
  return Sequence::SimData(rv.begin(),rv.end());
}
