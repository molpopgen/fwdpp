# Tutorial 3: Data serialization

## Intro

This tutorial covers the following topics:

* Writing simulated populations to files, and reading them back in
* Copying simulated population in memory

## The basic functions

For single-deme simulations, the library provides the following functions:

* fwdpp::write_binary_pop 
* fwdpp::read_binary_pop

These functions write and read simulated populations in a binary (_e.g._ not human-readable) format.

For multi-deme simulations, we have:

* fwdpp::write_binary_metapop 
* fwdpp::read_binary_metapop

As with the rest of the library, these functions are implemented using a combination of templates + function overloading.

Currently, these functions support the following types of simulations:

* single-locus, single-deme
* single-locus, multi-deme
* multi-locus, single-deme

The type requirements for these functions are:

* Containers of objects used in fwdpp-based simulations (containers of gametes, mutations, diploids, etc.)
* Output streams compatible with the public interface of [std::ostream](http://www.cplusplus.com/reference/ostream/ostream/)
* Input streams compatible with the public interfact of [std::istream](http://www.cplusplus.com/reference/istream/istream/)
* Input via the gzFile type defined by [zlib](http://zlib.net)

## Writing and reading mutations

You need to tell these functions how to read/write mutation objects.  Specifically, you need to:

* Define a mutation write function, which takes a mutation object and an output stream type as an object.
* Define a mutation read function, taking an input stream type as an object

For a concrete example of implementing "mutation readers" and "mutation writers", see the implementation of
fwdpp::mutation_writer and fwdpp::mutation_writer.

Below, you will see how to pass these to fwdpp::write_binary_pop and fwdpp::read_binary_pop.

## Reading and writing diploids.

Version 0.3.1 allowed the use of custom diploid genotype objects (@ref md_md_customdip).  Such objects can be serialized with no additional effort. _However, doing so will result in custom data associated with these custom types not being serialized._

In order to solve this, define your own custom structures for writing and reading:

~~~{.cpp}
struct diploid_writer {
	using result_type = void;
	template<typename dip_itr, typename streamtype >
	inline result_type( dip_itr i, streamtype & o ) const
	{
		//Do the right thing here.
	}
};
~~~

The arguments should be an iterator pointing to one of your custom diploids, and a reference to an output stream (including a gzFile!).  A function object of the same form should also be written to read the data back in.

These custom serialization objects may be passed as the last arguments to the functions discussed in the next section.  By default, they are implicitly passed fwdpp::diploidIOplaceholder, which is essentially an empty object.

## In-memory copying

It may be desirable to be able to restore a population from a previous state.  For example, you may wish to repeatedly introduce a mutation at a position, and then simulate until it is fixed.  For replicates where it is lost, you will restore the population to the state it was in when you introduced the mutation.  You may do this using fwdpp::write_binary_pop  or fwdpp::write_binary_metapop, and write the population to an in-memory buffer such as [std::ostringstream](http://www.cplusplus.com/reference/sstream/ostringstream/).  You may restore the population by reading the data from a [std::istringstream](http://www.cplusplus.com/reference/sstream/istringstream/).

Here is some pseudocode for an single-deme simulation:

~~~{.cpp}
//run a simulation...

//write it to a buffer:
std::ostringstream buffer;
fwdpp::write_binary_pop(gametes,mutations,diploids,std::bind(mwriter(),std::placeholders::_1,std::placeholders::_2),buffer);

//Now, if you want to restore it:
std::istringstream inbuffer(buffer.str());
decltype(mutations) mutations2;
decltype(gametes) gametes2;
decltype(diploids) diploids2;
fwdpp::read_binary_pop(gametes2,
			 mutations2,
			 diploids2,
			 std::bind(mreader(),std::placeholders::_1),
			 inbuffer);
~~~

This approach works because fwdpp::write_binary_pop and fwdpp::write_binary_metapop perform "deep copies" of the data, allowing the complete restoration of all the pointers, etc., stored in the containers.

## Writing to files

You may use the same functions to write/read to from files.  You may use and sort of [std::ostream](http://www.cplusplus.com/reference/ostream/ostream/)-compatible object for output.  You may use either a [std::istream](http://www.cplusplus.com/reference/istream/istream/)-compatible object for input or a gzFile object from [zlib](http://zlib.net).

## Managing output from multiple independent processes

In practice, we often implement simulations with the idea that one replicate will be run per process, and that we will spread a large number of processes across our cluster.  Typically, we collect the output of these simulations into a single file.  To do this, we implement POSIX-style advisory file locking.  See [here](https://github.com/molpopgen/BigDataFormats) for a discussion of how to do that using the low-level C functions found in <fcntl.h>.   A more modern C++-based approach would be to use the [boost](http://www.boost.org) synchronization library.  See [here](https://gist.github.com/molpopgen/651e4ac81253f34364f7) for simple examples.

### Why not use boost serialization?

Cliff's notes version: this approach would work well for a single process writing lots of replicates to one file.  However, the archives that this library creates do not support appending on to an existing library.  This limitation is unfortunate, as being able to use this library would simplify some of the implementation.
