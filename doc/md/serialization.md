# Tutorial 3: Data serialization

## Intro

This tutorial covers the following topics:

* Writing simulated populations to files, and reading them back in
* Copying simulated population in memory

## The basic functions

For single-deme simulations, the library provides the following functions:

* KTfwd::write_binary_pop 
* KTfwd::read_binary_pop

These functions write and read simulated populations in a binary (_e.g._ not human-readable) format.

For multi-deme simulations, we have:

* KTfwd::write_binary_metapop 
* KTfwd::read_binary_metapop

As with the rest of the library, these functions are implemented using a combination of templates + function overloading.

Currently, these functions support the following types of simulations:

* single-locus, single-deme, gamete-based 
* single-locus, multi-deme, gamete-based
* single-locus, single-deme, individual-based 
* single-locus, multi-deme, individual-based
* multi-locus, single-deme, individual-based 

The type requirements for these functions are:

* Containers of objects used in fwdpp-based simulations (containers of gametes, mutations, diploids, etc.)
* Output streams compatible with the public interface of [std::ostream](http://www.cplusplus.com/reference/ostream/ostream/)
* Input streams compatible with the public interfact of [std::istream](http://www.cplusplus.com/reference/istream/istream/)
* Input via the gzFile type defined by [zlib](http://zlib.net)

## Writing and reading mutations

You need to tell these functions how to read/write mutation objects.  Specifically, you need to:

* Define a mutation write function, which takes a mutation object and an output stream type as an object.
* Define a mutation read function, taking an input stream type as an object

Let's look at an example.  For this definition of a mutation type:

~~~{.cpp}
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
~~~

This function object works as a writer:

~~~{.cpp}
//function object to write mutation data in binary format
struct mwriter
{
  typedef void result_type;
  result_type operator()( const mutation_with_age & m, std::ostream & buffer ) const
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
~~~

And this function object works as a reader:

~~~{.cpp}
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
~~~

In the next section, you'll see how to pass these to KTfwd::write_binary_pop and KTfwd::read_binary_pop.

## In-memory copying

It may be desirable to be able to restore a population from a previous state.  For example, you may wish to repeatedly introduce a mutation at a position, and then simulate until it is fixed.  For replicates where it is lost, you will restore the population to the state it was in when you introduced the mutation.  You may do this using KTfwd::write_binary_pop  or KTfwd::write_binary_metapop, and write the population to an in-memory buffer such as [std::ostringstream](http://www.cplusplus.com/reference/sstream/ostringstream/).  You may restore the population by reading the data from a [std::istringstream](http://www.cplusplus.com/reference/sstream/istringstream/).

Here is some pseudocode for an individual-based simulation:

~~~{.cpp}
//run a simulation...

//write it to a buffer:
std::ostringstream buffer;
KTfwd::write_binary_pop(&gametes,&mutations,&diploids,std::bind(mwriter(),std::placeholders::_1,std::placeholders::_2),buffer);

//Now, if you want to restore it:
std::istringstream inbuffer(buffer.str());
decltype(mutations) mutations2;
decltype(gametes) gametes2;
decltype(diploids) diploids2;
KTfwd::read_binary_pop(&gametes2,
			 &mutations2,
			 &diploids2,
			 std::bind(mreader(),std::placeholders::_1),
			 inbuffer);
~~~

This approach works because KTfwd::write_binary_pop and KTfwd::write_binary_metapop perform "deep copies" of the data, allowing the complete restoration of all the pointers, etc., stored in the containers.

## Writing to files

You may use the same functions to write/read to from files.  You may use and sort of [std::ostream](http://www.cplusplus.com/reference/ostream/ostream/)-compatible object for output.  You may use either a [std::istream](http://www.cplusplus.com/reference/istream/istream/)-compatible object for input or a gzFile object from [zlib](http://zlib.net).

## Managing output from multiple independent processes

In practice, we often implement simulations with the idea that one replicate will be run per process, and that we will spread a large number of processes across our cluster.  Typically, we collect the output of these simulations into a single file.  To do this, we implement POSIX-style advisory file locking.  See [here](https://github.com/molpopgen/BigDataFormats) for a discussion of how to do that using the low-level C functions found in <fcntl.h>.   A more modern C++-based approach would be to use the [boost](http://www.boost.org) synchronization library.  See [here](https://gist.github.com/molpopgen/651e4ac81253f34364f7) for simple examples.

### Why not use boost serialization?

Cliff's notes version: this approach would work well for a single process writing lots of replicates to one file.  However, the archives that this library creates do not support appending on to an existing library.  This limitation is unfortunate, as being able to use this library would simplify some of the implementation.
