# Advanced topics: custom diploid genotypes

I should read this before going apeshit: http://stackoverflow.com/questions/14600201/why-should-i-avoid-stdenable-if-in-function-signatures

## Intro

This section discusses how to implement simulations where a diploid's genotype is represented by a user-defined type.

This is a document covering advanced programming issues using __fwdpp__.  Please see the introductory tutorials if you are new to programming with the library.

## Rationale

In all versions of __fwdpp__ \f$\leq\f$ 0.3.0, the genotype of a diploid (at a single locus, including a specific locus in a multilocus simulation) had the following type:

~~~{.cpp}
//Note: standard_diploid_t is not, nor ever was, a fwdpp type--this is simply a typedef for illustration's sake.
using standard_diploid_t = std::pair< glist::iterator, glist::iterator>;
~~~

And a single population was implemented as:

~~~{.cpp}
using standard_diploid_vec_t = std::vector< standard_diploid_t >;
~~~

Further, all fitness functions in the library (see KTfwd::site_dependent_fitness, for example) had the following form, and custom fitness models had to mimic this form as well (see @ref md_md_policies):

~~~{.cpp}
struct fitness_model {
	using result_type = double;
	template<typename iterator_t>
	inline operator()( const iterator_t & g1, const iterator_t & g2 ) const {
		//Do something useful and return a result_type
	}
};
~~~

This design meant that all diploids were _unlabelled_, meaning that things like separate sexes, geographic locations, etc., could not be modeled without writing explicit overloads of KTfwd::sample_diploid.  However, it _should_ be possible to implement some forms of these more complex simulations without having to overload the sampling functions.  Further, the library's template-based design should allow for these complexities to be dealt with at _compile-time_ instead of at run time.

In __fwdpp__ 0.3.1, I introduced the ability to specify diploid genotype types other than pairs of iterators to gametes.  The library maintains source code compatibility with existing simulations via a tag/dispatch model.  Basically, custom diploid genotype types are "tagged" at compile time.

## The minimal custom diploid type

Your custom diploid type is constrained in the following ways:

* Its gamete iterator variables must be called "first" and "second".
* The typedef first_type must exist, and be an alias for the iterator type.

These requirements are minimal, and force your custom type to have the same names/typedefs as the standard pair template.

Thus, the minimal valid custom diploid type (at least as far as GCC and clang++ are concerned) is:

~~~{.cpp}
//Your type must inherit from public KTfwd::tags::custom_diploid_t
struct diploid_t : public KTfwd::tags::custom_diploid_t
{
	using first_type = glist_t::iterator;
	//You don't need second_type here,
	//but I'm including it for completeness vis-a-vis std::pair
	using second_type = glist_t::iterator;
	//Iterator to gamete 1
	first_type first;
	//Iterator to gamete 2
	second_type second;
	//Constructor
	diploid_t() : first(first_type()),second(second_type()) {}
        diploid_t(first_type g1, first_type g2) : first(g1),second(g2){}
};
~~~

The base class is the "dispatch tag" alluded to in the previous section.  The "glist_t" is an alias for the doubly-linked list of gametes (see @ref md_md_datatypes).

### Separate sexes

We are now able to add more data to our diploids:

~~~{.cpp}
//Your type must inherit from public KTfwd::tags::custom_diploid_t
struct diploid_t : public KTfwd::tags::custom_diploid_t
{
	using first_type = glist_t::iterator;
	//You don't need second_type here,
	//but I'm including it for completeness vis-a-vis std::pair
	using second_type = glist_t::iterator;
	//Iterator to gamete 1
	first_type first;
	//Iterator to gamete 2
	second_type second;
	//"Sex" of individual
	bool female;
	//Constructor -- and you may want to write others in this case...
	diploid_t() : first(first_type()),second(second_type()),female(true) {}
        diploid_t(first_type g1, first_type g2) : first(g1),second(g2),female(true){}
};
~~~

## Fitness models

We may now consider fitness policies that depend on the additional data in our custom diploid types.  The fitness policies expecting two iterators to gametes are no longer sufficient, as they cannot know about the additional data in your custom type.  Thus, a fitness policy that depends on a custom diploid type must take an iterator pointing to the diploid genotype type as an argument.  In other words, a fitness policy must have the following form:

~~~{.cpp}
//Here, we assume diploid_t inherits from
// KTfwd::tags::custom_diploid_t
using diploid_vec_t = std::vector< diploid_t >;

struct my_new_fitness_pol {
	using result_type = double;
	inline result_type operator()( const diploid_vec_t::const_iterator & dip_itr ) const {
		//Do something interesting and return a double.
	}
};
~~~

In other words, the _signature_ of a fitness policy must be equivalent to:

~~~{.cpp}
std::function<double(const diploid_vec_t::const_iterator &)> my_new_fitness_pol = //something
~~~

It is straightforward to implement such policies as templates, too.  Let's look at a first pass at this:

~~~{.cpp}
struct my_new_fitness_pol_template {
	using result_type = double;
	template< typename diploid_genotype_itr >
	inline result_type operator()( const diploid_genotype_itr & dip_itr ) const {
	//Do something interesting and return a double.
	//dip_itr->first will access the first gamete, etc.
	}
};
~~~

The following built-in functions support custom diploids in this way:

* KTfwd::site_dependent_fitness
* KTfwd::additive_diploid
* KTfwd::multiplicative_diploid

### Fitness models in multilocus simulations

For multilocus simulations, a diploid is a vector of diploid genotype types, which by default is assumed to be a vector of pairs of iterators to gametes.  Thus, for multilocus sims involving custom types, there is no need to make use of the dispatch tags discussed above, unless you want to use custom diploid types _and_ the built-in fitness policies:

~~~{.cpp}
struct mloc_fitness {
  typedef double result_type;
  inline double operator()( const multiloc_t::dipvector_t::const_iterator & diploid ) const {
     using itr_t = multiloc_t::dipvector_t::const_iterator;
     //Fitness is additive across loci, and loci are additive over mutations:
     return std::accumulate( diploid.begin(), diploid.end(), 0, [](const double & w, const itr_t::value_type & __l ) {
     	    return w + KTfwd::additive_diploid(__l,2.);
	    } );
  }
};
~~~

## Cost vs. benefit

The pros of defining your own diploid type are:

* Being able to attach data to a diploid genotype beyond the two gametes that it contains.  One can imagine all sorts of things here, but the existing code base may not (yet) support them all.

The cons are:

* Fitness policies are written differently
* In order to use the sugar layer's streamlined methods for declaring population containers, you must pass your diploid type to those templates (see @ref md_md_sugar).

The cons mean that it is not trivial to switch an impementation back and forth between custom and non-custom diploid genotype representations.
