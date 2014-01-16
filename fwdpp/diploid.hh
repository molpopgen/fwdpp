/*! \file diploid.hh
  \brief Main header for programming using this library
  \warning Do not try to include individual headers a la carte. Life will get confusing for you. Just include this file.

  \code
  #include <fwdpp/diploid.hh>
  \endcode
 */
#ifndef __DIPLOID_HH__
#define __DIPLOID_HH__

//namespace std
#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <cassert>
#include <iterator>
#include <fstream>
#include <ctime>
#include <cmath>

//gsl
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

//boost
#include <boost/bind.hpp>

//headers from this project
#include <fwdpp/debug.hpp>
#include <fwdpp/migration.hpp>
#include <fwdpp/diploid_functions.hpp>
#include <fwdpp/fitness_models.hpp>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/fwd_functional.hpp>
#include <fwdpp/insertion_policies.hpp>
#include <fwdpp/sampling_functions.hpp>
#include <fwdpp/mutation.hpp>
#include <fwdpp/util.hpp>
#include <fwdpp/initms.hpp>
#include <fwdpp/IO.hpp>
#endif

/*! \mainpage Introduction
  \section purpose Purpose and Intent
  This page is the documentation for fwdpp, A C++ template library for efficient simulation of large populations with mutation, recombination, selection and migration.

  The library exists in namespace KTfwd.

  \subsection dependencies Dependencies On Other Libraries
  fwdpp depends on the following: \n
  1.  libsequence, available from www.molpopgen.org \n
  2.  The GNU Scientific library, available from www.gnu.org/software/gsl \n
  3.  The Boost C++ libraries, available from www.boost.org. \n

  \subsection examples How to use the library
  The best way to start will be to study the examples provided with the library and the \ref page1.  Regarding the examples, I would suggest doing so in the following order: \n
  1. diploid.cc: the most basic simulation of a Wright-Fisher population.  This is the most heavily-commented example.  Other examples assume that you've studied this one and understand it. \n
  2. diploid_binaryIO.cc: The same as diploid.cc but the output is changes to write the entire population to a file.  Mutiple processes run on a grid system can write to the same binary file due to use of POSIX file locking routines.\n
  3. diploid_fixed_sh.cc: The same as diploid.cc with the addition of a mutation rate to mutations with selection coefficient s and dominance h, with evolution under the multiplicative fitness model\n
  4. diploid_twopop_mig.cc: Evolution of a single population which then splits into two, followed by a period of migration.\n
  \n
  To use the library in your own code, include the library's main header:\n
  \code
  #include <fwdpp/diploid.hh>
  \endcode

  \subsection A note about compiling
  Several of the library functions allocate containers.  Based on performance testing that I have done, in some cases I use containers from namespace std, and in other cases I use the version found in namespace boost::container.  So, there is a hodge podge of container type usage within the library's internal code, and that is intentional, and the way it is set up results in the fastest run times on my system.   You have the option of reverting the code to using only the containers from namespace std by adding -DUSE_STANDARD_CONTAINERS when compiling programs using the library.  Note that the default use (e.g., leaving USE_STANDARD_CONTAINERS undefined) requires linking to the boost_system library (via -lboost_system).
 */

/*! \namespace KTfwd
  \brief The main (only) namespace defined by this library.
 */

/*! \page page1 Tutorial
  \section purpose Purpose
  This document will guide you through writing custom policies affecting the biology of a population in the simulation.

  \section model The model
  We will model the evolution of a locus where mutation effects are drawn from a Gaussian distribution with mean 0 and standard deviation sigma_s.\n
  \n
  The effect of a haplotype is the sum of mutation effects on that gamete, and the phenotypic value of a diploid is P = min(sum1,sum2) + R, where the two sums
  are the effects of haplotypes 1 and 2 and R is another Gaussian with mean 0 and standard deviation sigma_e.  The fitness of a diploid is then determined
  using a unit Gaussian with mean X. 
  \n
  \n
  Mutations will occur uniformly along the interval [0,1) according to the infinitely-many sites model, and recombination is uniform on the same interval.
  \n
  \section mutmodel Defining the mutation model
  There will be two mutation rates, u and us.  The former is the mutation rate per gamete to neutral mutations.  The latter is the mutation rate to selected mutations.
  Further, the position of a new mutation must not be the same as any currently-segregating mutation.\n
  \n
  We will make use of the existing mutation type KTfwd::mutation, which has a selection coefficient, s, added to the base class KTfwd::mutation_base.  For this model,
  haplotypes are the unit of evolution, so we can ignore the h (dominance) that is part of KTfwd::mutation, and we'll just leave it at the default value of 0.5.
  \n
  We will use a lookup table of current mutation positions in order to make sure that we assign novel positions to new mutations.  The most effective way to do that
  is to use a hash table.  I recommend boost::unordered_set.  It is used as follows:
  \code
  #include <boost/unordered_set.hpp>

  typedef boost::unordered_set<double,boost::hash<double>,KTfwd::equal_eps > lookup_table_type;
  \endcode

  In the above code block, we have replace the default std::operator== with KTfwd::equal_eps.  This will result in two doubles being considered equal if the
  difference between them is <= std::numeric_limits<double>::epsilon(), helping us skirt numeric precision issues inherent in the comparison of floating point values.\n
  \n
  Our mutation function will use the GNU Scientific Library for random number generation, and is implemented as follows:
  \code
  KTfwd::mutation custom_mutation_model( gsl_rng * r, const double & u, const double & us,const double & sigma_s,lookup_table_type * lookup)
  {
  double pos = gsl_rng_uniform(r);
  while( lookup->find(pos) != lookup->end() )
  {
  pos = gsl_rng_uniform(r);
  }
  lookup->insert(pos); //update the lookup table w/the accepted position
  double s = 0.;
  //ask if mutation is selected
  if( gsl_rng_uniform(r) <= us/(us+u) )
  {
  s = gsl_ran_gaussian(r,sigma_s);
  }
  //return a new mutation at position pos with effect s and count of 1 in the population.
  //Internally, the mutation constructor will check if s==0 or not set KTfwd::mutation::neutral to true or false as appropriate
  return KTfwd::mutation(pos,s,1);
  }
  \endcode
  In order to pass this mutation policy to KTfwd::mutate, we must construct a function call that KTfwd::mutate may access.  This is most easily accomplished via
  boost::bind:
  \code
  boost::bind( custom_mutation_model(r,u,us,sigma_s,&lookup) );
  \endcode
  \n
  As an aside, I typically pass a non-const pointer to the mutation list.  This allows me to do some extra checking to make sure that the lookup
  table is doing its job.  See diploid.cc for and example.

  \section fitnessmodel Defining the fitness model
  A fitness model will be passed iterators two gametes by one of the KTfwd::sample_diploid functions.  The programmer may define the fitness function to take
  other arguments and pass placeholders to the two gametes using boost::bind.  I typically defined fitness models as function objects rather than functions, which
  is convienient if they will be templates (because operator() can deduce template types during instantiation more readily than regular functions):

  \code
  struct custom_fitness_model
  //Yes, this may be a weird fitness model.  Dominance of a haplotype depends on the value of X, etc.
  {
  typedef double result_type;
  template<typename iterator>
  inline result_type operator()( const iterator & gamete1, const iterator & gamete2, gsl_rng * r, const double & sigma_e, const double & X ) const
  {
  result_type sum1=0.,sum2=0;
  for( unsigned i = 0 ; i < gametes1.smutations.size() ; ++i )
  {
  sum1 += i->s;
  }
  for( unsigned i = 0 ; i < gametes2.smutations.size() ; ++i )
  {
  sum2 += i->s;
  }
  result_type pheno = std::min(sum1,sum2) + gsl_ran_gaussian(r,sigma_e);
  //fitness is Gaussian with mean X and std_dev = 1
  result_type fitness = std::exp( (-1.*std::pow( (X-pheno),2.))/2. );
  return fitness;
  }
  };
  \endcode 
  To pass the fitness function to KTfwd::sample_diploid, use boost::bind:
  \code
  //_1 and _2 are placeholders for the iterators to gametes.  They will be passed in via sample_diploid.
  boost::bind( custom_fitness_model(_1,_2,r,sigma_e,X) );
  \endcode

  Please note that some C++ compilers may give errors with std::exp and std::pow even when <cmath> is included
  instead of <math.h>.  In such cases, replace with exp and pow, respectively.  Your mileage may vary.  Not
  all compilers are on the same page with putting the C function in namespace std.

  \section model2 Selection at sites linked to the sampled region.
  Here, we implement a mutation model for the recurrent hitch-hiking (RHH) model.  In this model, beneficial mutations
  enter the population at a constant rate.  These beneficial mutations have a constant selection coefficient and are 
  codominant with respect to fitness. We model the evolution of a region of size theta = 4Nu and rho=4Nr in which
  neutral and selected mutations occur.  We allow selected mutations to arise at positions linked to the region up
  to a genetic distance of s/rbp genetic units away (where rbp here is the recombination rate per base pair), following the classic theory.
  \n
  \n
  The relevant parameters of this model are:
  Lambda = the number of fixations of beneficial mutations per 4N generations per region of site theta
  r = the recombination rate per diploid per generation per site
  L = the length of the sampled region in base pairs.
  s = the selection coefficient 
  \n
  \n
  From the above parameters, and from classical population genetic theory, the mutation rate (per gamete, per site)
  to mutations with selection coefficient s is:\n
  \code
  //use the result on fixation probability for codominant mutation with selection coefficient s to work backwards to mutation rate per size
  double mu_pos_bp = (lambda > 0.) ? ((lambda*(1-exp(-4.*double(N0)*s)))/(8.*pow(double(N0),2.)*(1.-exp(-2.*s))))/double(nsites) : 0.;
  \endcode

  From theory, the maximum genetic distance at which a selected site will generate a hitch-hiking effect is:
  \code
  double maxd = 4*N*s;
  \endcode

  We will use maxd in the following way:  neutral mutations occur with position on the interval [0,1).  Selected sites will occur
  with positions uniformly-distributed between -maxd and 1+maxd, to allow them to occur both within the neutral region and up to a distance maxd
  on eiter side.  Further, the position of recombination events will be uniform on the interval [-maxd,1+maxd).

  The recombination rates are:
  \code
  //this is the recombination rate in the region of size theta
  double littler_diploid_region = rho/double(4*N);
  //this is the total recombination rate out to +/- maxd genetic units away from the region
  double littler_diploid_ttl = littler_diploid_region + 2.*maxd/rho
  \endcode
  We have the following other mutation rates:
  \code
  double mu_neutral = SOMETHING; //defined somehow, applies only to neutral region of size theta
  //selected mutation rate is selected mutation rate in region of size theta + selected mutation rate at flanking sites
  double selected_mu_ttl =mu_pos_bp*double(L) + mu_pos_bp*2.*maxd/(littler/double(L-1))
  //the probabilty that a selected mutation is within the sampled region is:
  double pwithin = mu_pos_bp*double(L)/selected_mu_ttl;
  \endcode

  We may wish to keep track of when mutations arise in the population, so that we can track fixations times, etc.:
  /code
  struct mutation_with_origin : public KTfwd::mutation
  //records the generation when the mutation entered the population
  {
  //o is the generation when the mutation arose (aka, the "origination time")
  mutable unsigned o;
  mutation_with_origin (const unsigned & __o ,const double & position, const double & sel_coeff,const unsigned & count,
  const double & dominance = 0.5) :
  o(__o),KTfwd::mutation(position,sel_coeff,count,dominance)
  {
  }
  };
  typedef mutation_with_origin mtype;
  \endcode
  We can now write our mutation model:
  \code
  typedef boost::unordered_set<double,boost::hash<double>,KTfwd::equal_eps > lookup_table_type;
  mtype RHH_mutation_model( gsl_rng * r, const unsigned & generation, const double mu_neutral, const double & mu_selected,
			  const double & p_selected_within,
			  const double & s, const double & h, const double & maxd, lookup_table_type * lookup ) //keep h (dominance) to 1 to match theory
{
  //Will the mutation be neutral, or a selected one?
  double neutral = (gsl_rng_uniform(r) <= mu_neutral/(mu_neutral+mu_selected)) ? true : false;
  double pos = std::numeric_limits<double>::max(); //set to some impossible dummy value
  double smut=0.;
  if( neutral )
    {
      pos = gsl_rng_uniform(r);
      while( lookup->find(pos) != lookup->end() )
	{
	  pos = gsl_rng_uniform(r);
	}
    }
  else
    {
      bool within = (gsl_rng_uniform(r) <= p_selected_within) ? true : false;
      pos = (within) ? gsl_rng_uniform(r) : gsl_ran_flat(r,0.,maxd);
      if( ! within )
	{
	  if (gsl_rng_uniform(r) <= 0.5 )
	    {
	      pos *= -1.;
	    }
	  else
	    {
	      pos += 1.;
	    }
	}
      while( lookup->find(pos) != lookup->end() )
	{
	  pos = (within) ? gsl_rng_uniform(r) : gsl_ran_flat(r,0.,maxd);
	  if( ! within )
	    {
	      if (gsl_rng_uniform(r) <= 0.5)
		{
		  pos *= -1.;
		}
	      else
		{
		  pos += 1.;
		}
	    }
	}
      smut=s;
    }
  lookup->insert(pos);
  return mtype(generation,pos,smut,1,h);
}
  \endcode

  And we can pass that to KTfwd::mutate using boost::bind:
  \code
  boost::bind(RHH_mutation_model,r,mu_neutral,mu_selected_ttl,pwithin,s,h,maxd,&lookup);
  \endcode

  Note that the above mutation model specifies a region from position from -maxd to maxd+1.
  The interval from 0 to 1 is where neutral and selected mutations arise.  The flanking bits 
  out to +/- maxd or so are where selected mutations arise.  Thus, we also need recombination
  on the interval from -maxd to maxd+1.  We pass this function to KTfwd::recombine as the genetic map:
  \code
  double recurrent_sweep_genetic_map(gsl_rng * r, const double & littler_neut,
				   const double & ttl_rec_rate,
				   const double & maxd)
{
  double rdm = gsl_rng_uniform(r),pos;
  if( rdm <= littler_neut/ttl_rec_rate )
    {
      pos = gsl_rng_uniform(r);
    }
  else
    {
      pos = (gsl_rng_uniform(r) <= 0.5) ? gsl_ran_flat(r,-1.*maxd,0.) : gsl_ran_flat(r,1.,1.+maxd);
    }
  return pos;
}
  \endcode

  The above is used as follows:
\code
	      unsigned nrec = KTfwd::recombine(r, 
					       &gametes,
					       twoN, 
					       littler,
					       //genetic map uniform over total region, neutral + selected
					       boost::bind( recurrent_sweep_genetic_map,r,littler_neut,
							    littler,maxd ) );
\endcode
  in order to get uniform recombination over the interval.

  The appropriate fitness model would be struct KTfwd::multiplicative_diploid with a scaling of 2.0.

The fully-worked out code for this example is in RHH.cc, which also demonstrates how to print out fixation times, etc.
*/
