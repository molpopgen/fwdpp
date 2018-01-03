# Writing data to binary formats

It is often useful to be able to serialize a population into a binary format and re-construct a serialized population.
The namespace fwdpp::io provides the following functions for these operations:

* fwdpp::io::serialize_population
* fwdpp::io::deserialize_population

## Required specializations

**Note** This section is only required if you use functions in namespace fwdpp::io.

The above functions are templates that require certain specializations to be visible at compile time in order to work
properly.  These specializations are of function objects to serialize/deserialize mutations, gametes, and diploids.
Their names are:

* fwdpp::io::serialize_mutation and fwdpp::io::deserialize_mutation
* fwdpp::io::serialize_gamete and fwdpp::io::deserialize_gamete
* fwdpp::io::serialize_diploid and fwdpp::io::deserialize_diploid

The default definitions of these types have the following behaviors:

* fwdpp::io::serialize_mutation and fwdpp::io::deserialize_mutation throw an exception.  Therefore, you *must*
  specialize them for your mutation type.  The built-in mutation types already have such specializations, which you may
  study if you develop your own mutation type.  See fwdpp::io::serialize_mutation\<fwdpp::popgenmut\> and
  fwdpp::io::deserialize_mutation\<fwdpp::popgenmut\> for examples.
* The gamete serialization objects assume the gamete type is fwdpp::gamete.  If this is not true, code may fail to
  compile and you need a specialization.
* The diploid serialization code only deals with the `first` and `second` data fields.  If your diploid type requires
  additional meta-data, then you need a specialization. See fwdpp::io::serialize_diploid\<custom_diploid_testing_t\> and
  fwdpp::io::deserialize_diploid\<custom_diploid_testing_t\>, which are part of the fwdpp unit test suite.

## Low-level functions

The functions fwdpp::io::serialize_population and fwdpp::io::deserialize_population discussed above are implemented in
terms of the following low-level functions:

* fwdpp::io::write_mutations and fwdpp::io::read_mutations
* fwdpp::io::write_gametes and fwdpp::io::read_gametes
* fwdpp::io::write_diploids and fwdpp::io::read_diploids

You may find these functions useful for more nuanced approaches to serialization.

## Streams

The stream template types all assume interfaces that are modeled after
[std::basic_iostream](http://en.cppreference.com/w/cpp/io/basic_iostream).  The boost
[iostreams](http://www.boost.org/doc/libs/1_66_0/libs/iostreams/doc/index.html) library will satisfy this criterion.

## General guidelines

Working with streams is a C++ topic rather than a fwdpp topic.  However, keep the following in mind:

* The [sstream](http://en.cppreference.com/w/cpp/header/sstream) header may be used to serialize in-memory.
* The objects in [sstream](http://en.cppreference.com/w/cpp/header/sstream) return a `std::string` via the `str()`
  member function.  Further, a `std::string`'s `c_str()` member function returns the underlying C-like `char *`.
* Most compression libraries ([zlib](http://zlib.net), [bzip2](http://bzip2.org)) expect a `char *`.

Therefore, you may serialize in-memory and then write the results to a compressed file.

## Limitations

### Large populations

When simulating very large populations, it is conceivalbe that you could run out of RAM if you attempt to serialize
in-memory.  One current workaround is to use the low-level object writing types (fwdpp::io::serialize_mutation *et al.*)
to roll your own serialization code that flushes buffers more often.  If you take a look at the current serialization
code, you'll see that doing this would not be a difficult task.  Another workaround may be to simply serialize to an
uncompressed temporary file.

Future library versions may provide built-in functions
to facilitate serialization in "chunks".

