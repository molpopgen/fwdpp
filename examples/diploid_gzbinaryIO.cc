/*
  \include diploid_gzbinaryIO.cc

  Same as diploid.cc, but this version uses fwdpp/IO.hpp to write the population to a gzipped binary file.
  
  The population is then read back in and compared to what was written out.

  Main point here is to show how the write/read function objects for fwdpp/IO.hpp should be written.

  Also illustrates POSIX file locking via <fcntl.h>, which is super-useful on clusters.

  Example use that runs quickly: ./diploid_gzbinaryIO 1000 10 10 10000 1 index haps.bin.gz $RANDOM
*/

#include <fwdpp/diploid.hh>
#include <Sequence/SimData.hpp>
#include <numeric>
#include <functional>
#include <cassert>

#include <fcntl.h>

#include <boost/unordered_set.hpp>
#include <boost/container/list.hpp>
#include <boost/container/vector.hpp>
#include <boost/pool/pool_alloc.hpp>

struct mutation_with_age : public KTfwd::mutation_base
{
  double s,h;
  unsigned g; //generation in which mutation arose
  mutation_with_age(const unsigned & __o,const double & position, const unsigned & count, 
		    const double & __s, const double __h,
		    const bool & isneutral = true)
    : KTfwd::mutation_base(position,count,isneutral),s(__s),h(__h),g(__o)
  {	
  }
};

//function object to write mutation data in binary format
struct mwriter
{
  typedef void result_type;
  result_type operator()( const mutation_with_age & m, std::ostringstream & buffer ) const
  {
    unsigned u = m.n;
    buffer.write( reinterpret_cast< char * >(&u),sizeof(unsigned) );
    u = m.g;
    buffer.write( reinterpret_cast< char * >(&u),sizeof(unsigned) );
    bool b = m.neutral;
    buffer.write( reinterpret_cast< char * >(&b),sizeof(bool) );
    double d = m.pos;
    buffer.write( reinterpret_cast< char * >(&d),sizeof(double) );
    d = m.s;
    buffer.write( reinterpret_cast< char * >(&d),sizeof(double) );
    d = m.h;
    buffer.write( reinterpret_cast< char * >(&d),sizeof(double) );
  }
};


/*
  Function object to read mutation data in binary format.
  This is the difference from diploid_binaryIO.cc,
  as it is based on gzread from a gzFile rather than
  an istream.

  Note: zlib.h is already included via fwdpp,
  but there would be no harm in doing so again above.
*/
struct mreader
{
  typedef mutation_with_age result_type;
  result_type operator()( gzFile gzin ) const
  {
    unsigned n;
    gzread( gzin,&n,sizeof(unsigned) );
    unsigned g;
    gzread( gzin,&g,sizeof(unsigned) );
    bool neut;
    gzread( gzin,&neut,sizeof(bool) );
    double pos;
    gzread( gzin,&pos,sizeof(double) );
    double s;
    gzread( gzin,&s,sizeof(double) );
    double h;
    gzread( gzin,&h,sizeof(double) );
    return result_type(g,pos,n,s,h,neut);
  }
};


typedef mutation_with_age mtype;
typedef boost::pool_allocator<mtype> mut_allocator;
typedef boost::container::list<mtype,mut_allocator > mlist;
typedef KTfwd::gamete_base<mtype,mlist> gtype;
typedef boost::pool_allocator<gtype> gam_allocator;
typedef boost::container::vector<gtype,gam_allocator > gvector;

typedef boost::unordered_set<double,boost::hash<double>,KTfwd::equal_eps > lookup_table_type;

mutation_with_age neutral_mutations_inf_sites(gsl_rng * r,const unsigned & generation,mlist * mutations,
					      lookup_table_type * lookup)
{
  double pos = gsl_rng_uniform(r);
  while( lookup->find(pos) != lookup->end() )
    {
      pos = gsl_rng_uniform(r);
    }
  lookup->insert(pos);
  assert(std::find_if(mutations->begin(),mutations->end(),boost::bind(KTfwd::mutation_at_pos(),_1,pos)) == mutations->end());
  return mutation_with_age(generation,pos,1,0.,0.,true);
}

int main(int argc, char ** argv)
{
  if (argc != 9)
    {
      std::cerr << "Too few arguments\n"
		<< "Usage: diploid_gzbinaryIO N theta rho ngens replicate_no indexfile hapfile seed\n";
      exit(10);
    } 
  int argument=1;
  const unsigned N = atoi(argv[argument++]);
  const double theta = atof(argv[argument++]);
  const double rho = atof(argv[argument++]);
  const unsigned ngens = atoi(argv[argument++]);
  const unsigned replicate_no = atoi(argv[argument++]);
  const char * indexfile = argv[argument++];
  const char * hapfile = argv[argument++];
  const unsigned seed = atoi(argv[argument++]);

  const double mu = theta/double(4*N);
  const double littler = rho/double(8*N);
  
  std::copy(argv,argv+argc,std::ostream_iterator<char*>(std::cerr," "));
  std::cerr << '\n';
  
  gsl_rng * r =  gsl_rng_alloc(gsl_rng_ranlxs2);
  gsl_rng_set(r,seed);

  //the population begins with 1 gamete with no mutations

  unsigned twoN = 2*N;

  gvector gametes(1,gtype(twoN));
  mlist mutations;
  std::vector<mtype> fixations;
  std::vector<unsigned> fixation_times;
  unsigned generation;

  fixations.clear();
  fixation_times.clear();
  double wbar;

  lookup_table_type lookup;
  for( generation = 0; generation < ngens; ++generation )
    {
      wbar = KTfwd::sample_diploid(r,&gametes,twoN,
				   boost::bind(KTfwd::multiplicative_diploid(),_1,_2,2.),
				   boost::bind(KTfwd::mutation_remover(),_1,0,twoN));
      KTfwd::remove_fixed_lost(&mutations,&fixations,&fixation_times,&lookup,generation,twoN);
      assert(KTfwd::check_sum(gametes,twoN));
      assert( lookup.size() == mutations.size() );
      unsigned nmuts = KTfwd::mutate(r,&gametes,&mutations,mu,
				     boost::bind(neutral_mutations_inf_sites,r,generation,_1,&lookup),
				     boost::bind(KTfwd::push_at_end<gtype,gvector >,_1,_2),
				     boost::bind(KTfwd::insert_at_end<mtype,mlist>,_1,_2));

      assert(KTfwd::check_sum(gametes,twoN));
      unsigned nrec = KTfwd::recombine(r, &gametes, twoN, littler, boost::bind(gsl_rng_uniform,r));
      assert(KTfwd::check_sum(gametes,twoN));
    }
  std::ostringstream buffer;
      
  KTfwd::write_binary_pop(&gametes,&mutations,boost::bind(mwriter(),_1,_2),buffer);

  FILE * fp = fopen( hapfile, "ab" );
  fseek(fp, -4, SEEK_END);
  std::cerr << ftell(fp) << '\n';
  uint32_t gzsize;
  fread(&gzsize,4,1,fp);
  std::cerr << ftell(fp) << '\n';
  std::cerr << "size = " << gzsize << '\n';
  z_off_t offset = ftell(fp);
  gzFile gzf = gzopen( hapfile, "ab" ); //append mode, binary

  //We now need to get the size of the gzfile.
  //z_off_t offset = gzoffset(gzf);
  std::cerr << "offset = " << offset << '\n';
  gzwrite( gzf, buffer.str().c_str(), buffer.str().size() );
  std::cerr << gzoffset( gzf ) << ' ' << gztell( gzf ) << '\n';
  gzclose( gzf );


  //now, read the data back in...
  std::vector< gtype > gametes2;
  mlist mutations2;

  gzf = gzopen( hapfile, "rb" ); //read mode, binary
  gzseek( gzf, offset, SEEK_SET );
  read_binary_pop(&gametes2,&mutations2,boost::bind(mreader(),_1),gzf);
  gzclose(gzf);
  //Now, compare what we wrote to what we read
  std::cout << "Mutations:\n";
  mlist::iterator mitr1 = mutations.begin(),mitr2=mutations2.begin();
  for( ; mitr1 != mutations.end() ; ++mitr1,++mitr2 )
    {
      std::cout << mitr1->pos << '\t' << mitr2->pos << '\t'
		<< mitr1->n << '\t' << mitr2->n << '\n';
    }
  std::cout << "Gametes:\n";
  for( unsigned i = 0 ; i < gametes.size() ; ++i )
    {
      std::cout << i << ':' << gametes[i].n << ' ' << gametes[i].mutations.size() << ' ';
      for(unsigned j = 0 ; j < gametes[i].mutations.size() ; ++j)
	{
	  std::cout << '(' 
		    << gametes[i].mutations[j]->pos << ','
		    << gametes[i].mutations[j]->n << ','
		    << gametes[i].mutations[j]->neutral 
		    << ')';
	}
      std::cout << '\n';

      std::cout << i << ':' << gametes2[i].n << ' ' << gametes2[i].mutations.size() << ' ';
      for(unsigned j = 0 ; j < gametes2[i].mutations.size() ; ++j)
	{
	  std::cout << '(' 
		    << gametes2[i].mutations[j]->pos << ','
		    << gametes2[i].mutations[j]->n << ','
		    << gametes2[i].mutations[j]->neutral 
		    << ')';
	}
      std::cout << '\n';
    }
}
