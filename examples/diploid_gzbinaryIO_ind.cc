/*
  \include diploid_gzbinaryIO_ind.cc

  Same as diploid_binaryIO_ind.cc, but input/output are gzip compressed
  
  The population is then read back in and compared to what was written out.

  Main point here is to show how the write/read function objects for fwdpp/IO.hpp should be written.

  Also illustrates POSIX file locking via <fcntl.h>, which is super-useful on clusters.

  Example use that runs quickly: ./diploid_hzbinaryIO 1000 10 10 10000 1 index haps.bin.gz $RANDOM
*/

#include <fwdpp/diploid.hh>
#include <Sequence/SimData.hpp>
#include <numeric>
#include <functional>
#include <cassert>

#include <fcntl.h>

#include <boost/bind.hpp>
#include <boost/function.hpp>
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
  This is the difference from diploid_binaryIO_ind.cc,
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
typedef boost::container::list<gtype,gam_allocator > glist;

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
		<< "Usage: diploid_gzbinaryIO_ind N theta rho ngens replicate_no indexfile hapfile seed\n";
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
  const double littler = rho/double(4*N);
  
  std::copy(argv,argv+argc,std::ostream_iterator<char*>(std::cerr," "));
  std::cerr << '\n';
  
  gsl_rng * r =  gsl_rng_alloc(gsl_rng_ranlxs2);
  gsl_rng_set(r,seed);

  //the population begins with 1 gamete with no mutations

  unsigned twoN = 2*N;

  glist gametes(1,gtype(twoN));
  mlist mutations;
  std::vector< std::pair< glist::iterator,glist::iterator > > diploids(N,std::make_pair(gametes.begin(),
											gametes.begin()));
  std::vector<mtype> fixations;
  std::vector<unsigned> fixation_times;
  unsigned generation;

  fixations.clear();
  fixation_times.clear();
  double wbar;

  lookup_table_type lookup;  

  //recombination map is uniform[0,1)  
  boost::function<double(void)> recmap = boost::bind(gsl_rng_uniform,r);

  for( generation = 0; generation < ngens; ++generation )
    {
      wbar = KTfwd::sample_diploid(r,
				   &gametes, 
				   &diploids,
				   &mutations,
				   N,     
				   mu,   
				   boost::bind(neutral_mutations_inf_sites,r,generation,_1,&lookup),
				   boost::bind(KTfwd::genetics101(),_1,_2,
					       &gametes,
					       littler,
					       r,
					       recmap),
				   boost::bind(KTfwd::insert_at_end<mtype,mlist>,_1,_2),
				   boost::bind(KTfwd::insert_at_end<gtype,glist>,_1,_2),
				   boost::bind(KTfwd::multiplicative_diploid(),_1,_2,2.),
				   boost::bind(KTfwd::mutation_remover(),_1,0,2*N));
      KTfwd::remove_fixed_lost(&mutations,&fixations,&fixation_times,&lookup,generation,twoN);
    }
  std::ostringstream buffer;
      
  KTfwd::write_binary_pop(&gametes,&mutations,&diploids,boost::bind(mwriter(),_1,_2),buffer);

  /*
    Note: gzFiles are not as easily compatible with file-locking and creating an index, etc.,
    as done in diploid_binaryIO.cc.  See https://github.com/molpopgen/BigDataFormats
    for why not.
   */

  //establish POSIX file locks for output
  struct flock index_flock;
  index_flock.l_type = F_WRLCK;/*Write lock*/
  index_flock.l_whence = SEEK_SET;
  index_flock.l_start = 0;
  index_flock.l_len = 0;/*Lock whole file*/

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

  //OK, we no have an exclusive lock on the index file.  In principle, that is sufficient for us to move on.

  //Now, write output to a gzipped file
  gzFile gzout = gzopen( hapfile, "ab" ); //open in append mode.  The b = binary mode.  Not required on all systems, but never hurts.
  gzwrite( gzout, buffer.str().c_str(), buffer.str().size() ); //write buffer to gzfile
  //write info to index file
  fprintf( index_fh, "%u %lld\n",replicate_no,gztell( gzout ) );
  gzclose(gzout);

  //We can now release close the index file, release the lock, etc.
  index_flock.l_type = F_UNLCK;
  if (fcntl(index_fd, F_UNLCK,&index_flock) == -1) 
    {
      std::cerr << "ERROR: could not release lock on " << indexfile << '\n';
      exit(10);
    }
  fflush( index_fh );
  fclose(index_fh);

   /*
    Read the data back in.
    Unlike our other examples, we must keep track of the cumulative offsets in the index file
    up until we find our record_number
  */
  std::ifstream index_in(indexfile);
  unsigned long long rec_offset = 0,offset_i;
  unsigned rid;
  while(! index_in.eof() )
    {
      index_in >> rid >> offset_i >> std::ws;
      if ( rid != replicate_no )
	{
	  rec_offset += offset_i;
	}
      else
	{
	  break;
	}
    }
  index_in.close();
  if(rid != replicate_no)
    {
      std::cerr << "Error: replicate id not found in index file, " << replicate_no << ' ' << rid << '\n';
      exit(10);
    }

  mlist mutations2;
  glist gametes2;
  std::vector< std::pair< glist::iterator,glist::iterator > > diploids2;

  gzFile gzin = gzopen( hapfile, "rb" ); //open it for reading.  Again, binary mode.
  KTfwd::read_binary_pop(&gametes2,
			 &mutations2,
			 &diploids2,
			 boost::bind(mreader(),_1),
			 gzin);

  //Now, compare what we wrote to what we read
  std::cout << gametes.size() << ' ' << gametes2.size() << ' ' << mutations.size() << ' ' << mutations2.size() 
	    << ' ' << diploids.size() << ' ' << diploids2.size() << '\n';

  for( unsigned i = 0 ; i < diploids.size() ; ++i )
    {
      std::cout << "Diploid " << i << ":\nWritten:\n";
      for( unsigned j = 0 ; j < diploids[i].first->mutations.size() ; ++j )
	{
	  std::cout << '(' << diploids[i].first->mutations[j]->pos << ','
		    << diploids[i].first->mutations[j]->n << ')';
	}
      std::cout << '\n';
      for( unsigned j = 0 ; j < diploids[i].second->mutations.size() ; ++j )
	{
	  std::cout << '(' << diploids[i].second->mutations[j]->pos << ','
		    << diploids[i].second->mutations[j]->n << ')';
	}
      std::cout << '\n';

      std::cout << "Read:\n";
      for( unsigned j = 0 ; j < diploids2[i].first->mutations.size() ; ++j )
	{
	  std::cout << '(' << diploids2[i].first->mutations[j]->pos << ','
		    << diploids2[i].first->mutations[j]->n << ')';
	}
      std::cout << '\n';
      for( unsigned j = 0 ; j < diploids2[i].second->mutations.size() ; ++j )
	{
	  std::cout << '(' << diploids2[i].second->mutations[j]->pos << ','
		    << diploids2[i].second->mutations[j]->n << ')';
	}
      std::cout << '\n';
    }
}
