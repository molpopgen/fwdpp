# Tutorial 2: Implementing multilocus simulations

In order to simulate discontiguous genomic segments, one has two options within __fwdpp__:

1. Write mutation, recombination, and fitness policies that "do the right thing" for your model.  For example, a mutation policy would need to know about the mutation rate at each locus, and appropriately assign mutations with the correct positions, fitness effects, etc.
2. Write separate mutation and recombination policies for each locus, and a fitness policy that calculates the fitness of a diploid over all loci.

This document covers the latter method, which I call the "multilocus" part of __fwdpp__.    I won't give any examples of the former method, as I'm opposed to the idea of having to develop, debug, and mainting large complex policies.  But if you want to do things that way, read the [tutorial on policies](@ref md_md_policies), because what you are looking for is all in there.

# Simple policies using the multi-locus machinery

__fwdpp__ only supports the implementation of individual-based multilocus simulations, and it is unlikely that these features will be added to the gamete-based part of the library.

The main conceptual difference between this part of the library and the examples shown in the [tutorial on policies](@ref md_md_policies) is the following:

* Instead of a single mutation model policy, you implement one mutation model per "locus".  These policies are stored in a vector and passed to KTfwd::sample_diploid.
* Similarly, you implement a recombination policy per locus, and pass a vector of those policies along to KTfwd::sample_diploid.
* Instead of a single list of gametes, you have a vector of lists of gametes, where each list represents the current gametes at a particular locus.
* A diploid is now represented as a vector of pairs of iterators derived from the vector of lists of gametes.
* A fitness policy calculates individual fitnesses from that vector of pairs of iterators.

At this point, it may be most useful to look at a concrete example.  The program diploid_ind_2locus.cc is distributed with the library source code, and we'll break down its essential parts in the next few sections.

### A multilocus mutation model

For this example, we simulate only neutral mutations, and the model has two loci.  Our function to return a new mutation will do the following:

1. Take the start position of a locus as an argument.
2. Return a new mutations whose position is uniformly-distributed on in interval \f$(\mathrm{start},\mathrm{start}+1]\f$

This is our mutation class:

~~~{.cpp}
struct mutation_with_age : public KTfwd::mutation_base
{
  unsigned g;
  double s,h;
  /*
    The constructor initializes all class data, including that of the base class via a constructor
    call to the base class.
  */
  mutation_with_age(const unsigned & __o,const double & position, const unsigned & count, const bool & isneutral = true)
    : KTfwd::mutation_base(position,count,isneutral),g(__o),s(0.),h(0.)
  {	
  }
};
~~~

It contains info for selection coefficients, etc., but we won't be using any of that.  (It is only there b/c I use the same mutation object for all of the example programs...)

Given the requirements outlined above, a function to return a mutation at a specific locus is:

~~~{.cpp}
//"beg" is the start position of this locus
mutation_with_age neutral_mutations_inf_sites(gsl_rng * r,const unsigned * generation,mlist * mutations,
					      lookup_table_type * lookup, const double & beg)
{
  //Generate new mutation position on the interval [0,1)
  double pos = gsl_ran_flat(r,beg,beg+1.);
  while( lookup->find(pos) != lookup->end() ) //make sure it doesn't exist in the population
    { 
      pos = gsl_ran_flat(r,beg,beg+1.);
    }
  //update the lookup table
  lookup->insert(pos);

  //In absence of DEBUG, make sure lookup table is working
  assert(std::find_if(mutations->begin(),mutations->end(),std::bind(KTfwd::mutation_at_pos(),std::placeholders::_1,pos)) == mutations->end());

  //return constructor call to mutation type
  return mutation_with_age(*generation,pos,1,true);
}
~~~

The details of the lookup table are covered in the main [tutorial](@ref md_md_policies) on policies.  It serves to efficiently ensure that we are sampling new mutation positions from an infinitely-many sites model.

We may now create mutation model policies by synthesizing a function call that returns a mutation_with_age using [std::bind](http://en.cppreference.com/w/cpp/utility/functional/bind):

~~~{.cpp}
auto mmodel0 = std::bind(neutral_mutations_inf_sites,r,&generation,std::placeholders::_1,&lookup,0.);
auto mmodel1 = std::bind(neutral_mutations_inf_sites,r,&generation,std::placeholders::_1,&lookup,1.);
~~~

The placeholder is for a pointer to a doubly-linked list of mutations (see argument 3 in the definition of neutral_mutations_inf_sites above).  Each of these is identical to what we do in a single-locus simulation.  The only "trick" here is that we pass a 0 to mmodel0, and a 1 to mmodel1 to represent the starting position for each locus, thus guaranteeing that neutral_mutations_inf_sites returns a mutation with positions \f$(0,1]\f$ or \f$(1,2]\f$.

We need to store these policies in a vector.  The easiest way to do that is to use the C++11 features [decltype](http://en.cppreference.com/w/cpp/language/decltype) and [list initialization](http://en.cppreference.com/w/cpp/language/list_initialization):

~~~{.cpp}
std::vector<decltype(mmodel0)> mmodels { mmodel0, mmodel1 };
~~~

### Digression: type signatures for for policies

The multilocus API requires that the user pass a vector of policies.  A vector's interface further requires that all elements contained by a vector are of the same type.  In the previous section, we used C++11's [auto](http://en.cppreference.com/w/cpp/language/auto) keyword to force the compiler to figure out the type of mmodel0 and mmodel1.  It just so happens that they are of the same type, which I think should be obvious if we consider that both are calls to std::bind with the same number of arguments, all of which are the same type.

Some users may try to implement the mutation models using C++11 [lambda expressions](http://en.cppreference.com/w/cpp/language/lambda):

~~~{.cpp}
auto mmodel0 = [&]( mlist * m ) { return neutral_mutations_inf_sites(r,&generation,m,&lookup,0.); };
auto mmodel1 = [&]( mlist * m ) { return neutral_mutations_inf_sites(r,&generation,m,&lookup,1.); };
std::vector< decltype(mmodel0) > mmodels { mmodel0, mmodel1 };
~~~

However, the above approach will fail because C++11 lambda expressions always have different types, even if their signatures are the same, as in the case above (both lambdas take the same number of arguments, all arguments are the same type, and the return value is the same).  In order to use lambda expressions for our mutation models, our vector of policies will have to specify a particular function signature for its object type using [std::function](http://en.cppreference.com/w/cpp/utility/functional/function):

~~~{.cpp}
std::vector< std::function<mutation_with_age(mlist *)> > mmodels {mmodel0,mmodel1};
~~~

Personally, I prefer the std::bind/decltype idiom here over the lambda/std::function approach, because it means less work on my part.


### Separate within-locus recombination policies

### Recombination between loci

### A (trivial) mutilocus fitness model

This needs special mention

## Mechanics of the multilocus recombination

Document the logic here
