/*
  \include diploid_binaryIO.cc

  Same as diploid.cc, but this version uses fwdpp/IO.hpp to write the population to a binary file.
  
  The population is then read back in and compared to what was written out.

  Main point here is to show how the write/read function objects for fwdpp/IO.hpp should be written.

  Also illustrates POSIX file locking via <fcntl.h>, which is super-useful on clusters.

  Example use that runs quickly: ./diploid_binaryIO 1000 10 10 10000 1 index haps.bin $RANDOM
*/

#include <fwdpp/diploid.hh>
#include <Sequence/SimData.hpp>
#include <numeric>
#include <functional>
#include <cassert>
#include <sstream>
#include <fcntl.h>

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

typedef mutation_with_age mtype;
#include <common_gamete.hpp>


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

//function object to read mutation data in binary format
struct mreader
{
  typedef mutation_with_age result_type;
  result_type operator()( std::istream & in ) const
  {
    unsigned n;
    in.read( reinterpret_cast< char * >(&n),sizeof(unsigned) );
    unsigned g;
    in.read( reinterpret_cast< char * >(&g),sizeof(unsigned) );
    bool neut;
    in.read( reinterpret_cast< char * >(&neut),sizeof(bool) );
    double pos;
    in.read( reinterpret_cast< char * >(&pos),sizeof(double) );
    double s;
    in.read( reinterpret_cast< char * >(&s),sizeof(double) );
    double h;
    in.read( reinterpret_cast< char * >(&h),sizeof(double) );
    return result_type(g,pos,n,s,h,neut);
  }
};


mutation_with_age neutral_mutations_inf_sites(gsl_rng * r,const unsigned & generation,mlist * mutations,
					      lookup_table_type * lookup)
{
  double pos = gsl_rng_uniform(r);
  while( lookup->find(pos) != lookup->end() )
    {
      pos = gsl_rng_uniform(r);
    }
  lookup->insert(pos);
  assert(std::find_if(mutations->begin(),mutations->end(),std::bind(KTfwd::mutation_at_pos(),std::placeholders::_1,pos)) == mutations->end());
  return mutation_with_age(generation,pos,1,0.,0.,true);
}

int main(int argc, char ** argv)
{
  if (argc != 9)
    {
      std::cerr << "Too few arguments\n"
		<< "Usage: diploid_binaryIO N theta rho ngens replicate_no indexfile hapfile seed\n";
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
				   std::bind(KTfwd::multiplicative_diploid(),std::placeholders::_1,std::placeholders::_2,2.),
				   std::bind(KTfwd::mutation_remover(),std::placeholders::_1,0,twoN));
      KTfwd::remove_fixed_lost(&mutations,&fixations,&fixation_times,&lookup,generation,twoN);
      assert(KTfwd::check_sum(gametes,twoN));
      assert( lookup.size() == mutations.size() );
      unsigned nmuts = KTfwd::mutate(r,&gametes,&mutations,mu,
				     std::bind(neutral_mutations_inf_sites,r,generation,std::placeholders::_1,&lookup),
				     std::bind(KTfwd::push_at_end<gtype,gvector >,std::placeholders::_1,std::placeholders::_2),
				     std::bind(KTfwd::insert_at_end<mtype,mlist>,std::placeholders::_1,std::placeholders::_2));

      assert(KTfwd::check_sum(gametes,twoN));
      unsigned nrec = KTfwd::recombine(r, &gametes, twoN, littler, std::bind(gsl_rng_uniform,r));
      assert(KTfwd::check_sum(gametes,twoN));
    }
  std::ostringstream buffer;
      
  KTfwd::write_binary_pop(&gametes,&mutations,std::bind(mwriter(),std::placeholders::_1,std::placeholders::_2),buffer);

  //establish POSIX file locks for output
  struct flock index_flock, hapfile_flock;
  index_flock.l_type = F_WRLCK;/*Write lock*/
  index_flock.l_whence = SEEK_SET;
  index_flock.l_start = 0;
  index_flock.l_len = 0;/*Lock whole file*/
      
  //lock index file ASAP
  FILE * index_fh = fopen(indexfile,"a");
  int index_fd = fileno(index_fh);
  if ( index_fd == -1 ) 
    { 
      std::cerr << "ERROR: could not open " << indexfile << '\n';
      exit(10);
    }
  if (fcntl(index_fd, F_SETLKW,&index_flock) == -1) 
    {
      std::cerr << "ERROR: could not obtain lock on " << indexfile << '\n';
      exit(10);
    }
      
  hapfile_flock.l_type = F_WRLCK;/*Write lock*/
  hapfile_flock.l_whence = SEEK_SET;
  hapfile_flock.l_start = 0;
  hapfile_flock.l_len = 0;/*Lock whole file*/

  FILE * haps_fh = fopen(hapfile,"a");
  int hapfile_fd = fileno(haps_fh);
      
  if ( hapfile_fd == -1 ) 
    { 
      std::cerr << "ERROR: could not open " << hapfile << '\n';
      exit(10);
    }
  if (fcntl(index_fd, F_SETLKW,&index_flock) == -1) 
    {
      std::cerr << "ERROR: could not obtain lock on " << hapfile << '\n';
      exit(10);
    }

  long int offset = ftell(haps_fh);
  //write the index
  fprintf( index_fh,"%ud\t%ld\n",replicate_no,offset );

  //write the buffered haplotype data
  write( hapfile_fd, buffer.str().c_str(), buffer.str().size() );

  //release locks and close files
  index_flock.l_type = F_UNLCK;
  hapfile_flock.l_type = F_UNLCK;

  if (fcntl(hapfile_fd, F_UNLCK,&hapfile_flock) == -1) 
    {
      std::cerr << "ERROR: could not releaselock on " <<  hapfile << '\n';
      exit(10);
    }
  fflush( haps_fh );
  fclose(haps_fh);
      
  if (fcntl(index_fd, F_UNLCK,&index_flock) == -1) 
    {
      std::cerr << "ERROR: could not release lock on " << indexfile << '\n';
      exit(10);
    }
  fflush( index_fh );
  fclose(index_fh);

  //now, read the data back in...
  std::vector< gtype > gametes2;
  mlist mutations2;

  std::ifstream in(hapfile,std::ios_base::in|std::ios_base::binary);
  in.seekg(offset);
  read_binary_pop(&gametes2,&mutations2,std::bind(mreader(),std::placeholders::_1),in);

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
