# Implementing multilocus simulations

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

### Separate recombination policies per locus

### A (trivial) mutilocus fitness model

This needs special mention

## Mechanics of the multilocus recombination

Document the logic here
