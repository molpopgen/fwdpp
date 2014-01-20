#ifndef __RECOMBINATION_METHODS_IND_HPP__
#define __RECOMBINATION_METHODS_IND_HPP__

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

namespace KTfwd
{

struct genetics101
/*! Genetics 101:
  Recombination results in 2 parental and 2 non-parental gametes.  
  Each of those 4 gametes is equally-likely to be inherited.
  Thus, with prob. of 50%, the offspring inherits a non-parental from parent 1,
  and likewise from parent 2, again with prob. of 50%.
  
  This means we don't have to run the recombination function all of the time.
  Rather, we just ask if the descendant inherets a parental or non-parental type
  from each of the parents.
  
  Also, if a parent is homozygous (e.g., p1g1 == p1g2), then we don't
  need to run the recombination function as parental and non-parental gametes 
  are identical.
  
  The above is useful b/c recombination can be computationally expensive

  This routine is probably appropriate when the total recombination distance is on the smaller side

  Sadly, details of multiple crossovers in different species are lacking.

  From JRI:

  2-DIMENSIONAL SPREADS OF SYNAPTONEMAL COMPLEXES FROM SOLANACEOUS PLANTS .6. HIGH-RESOLUTION RECOMBINATION NODULE MAP FOR TOMATO (LYCOPERSICON-ESCULENTUM)
  By: SHERMAN, JD; STACK, SM
  GENETICS  Volume: 141   Issue: 2   Pages: 683-708 

Hi Jeff,
 
Your original question was whether a single CO is more likely to occur between any particular pair of chromatids.  Other than a significant preference for using non-sister chromatids to repair meiotic double-strand breaks (as opposed to sister-chromatids) I’m not aware of any literature that demonstrates a bias (I’d be interested to see any you’ve found). 
 
Wayne and Kelly did a nice job of expanding on your original question to discuss interference between two COs.  There is a rich literature about CO interference along chromosomes (I’ll shamelessly point you to a review written by my lab - PMID: 20885817).  People have also looked for interference between chromatids.  In other words, if you have a CO between 1&3 are you more or less likely to see another one using the same two chromatids.  The naïve assumption is that the distribution of double COs among chromatids will follow a 1:2:1 ratio of two-strand:three-strand:four-strand DCOs (strands in this unfortunate lingo refers to chromatids).  We actually looked at this in Arabidopsis (PMID: 9419361) and found no divergence from this null-hypothesis.  If there was a bias (or chromatid interference) you’d expect to see skew away from the expected ratio.
 
I hope that helps.
 
If you want to dive into the very deep end of this area you might also want to ask Frank Stahl.
 
-cheers
 
Greg
 
*/
{
  typedef unsigned result_type;
  template<typename gamete_iterator_type,
	   typename gamete_list_type_allocator,
	   template<typename,typename> class gamete_list_type,
	   typename rec_pos_generator>
  unsigned operator()( gamete_iterator_type & g1,
		       gamete_iterator_type & g2,
		       gamete_list_type< typename gamete_iterator_type::value_type, gamete_list_type_allocator > * gametes,
		       const double & littler,
		       gsl_rng * r,
		       const rec_pos_generator & rp) const
  {
    typedef gamete_list_type< typename gamete_iterator_type::value_type, gamete_list_type_allocator > glist_t;
    unsigned NREC = 0;
    if( g1 != g2 && gsl_rng_uniform(r) <= 0.5 )
      //then a non-parental type is inherited from p1 and p1 has two different gametes
      {
	NREC += recombine_gametes(r,littler,gametes,g1,g2,rp,
				  boost::bind(update_if_exists_insert<typename gamete_iterator_type::value_type,glist_t>,_1,gametes));	  
      }
    return NREC;
  }	   
};

}

#endif //__RECOMBINATION_METHODS_IND_HPP__
