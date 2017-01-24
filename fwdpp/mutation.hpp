#ifndef _FWDPP_MUTATION_HPP_
#define _FWDPP_MUTATION_HPP_

#include <fwdpp/forward_types.hpp>
#include <fwdpp/fwd_functional.hpp>
#include <gsl/gsl_rng.h>

namespace KTfwd
{
    /*! \brief Apply mutation model to an individual gamete.  Used for
      individual-based forward simulations
      \param recycling_bin Recycling bin for mutations
      \param gamete_recycling_bin Recycling bin for gametes
      \param r GSL random number generator
      \param gametes Container of gametes
      \param mutations Container of mutations
      \param g gametes[g] is the gamete to be mutated
      \param mu the TOTAL mutation rate per gamete
      \param mmodel Mutation model policy
      \param gpolicy Policy determining how new gametes are added to population
      \param mpolicy Policy determining how new mutations are added to the
      population

      \note g is passed non-const and will be modified by mutation events.
      \note Used in invididual-based forward simulations.
      \return The location of the newly-mutated gamete in gametes.
    */
    template <typename queue_type, typename queue_type2,
              typename mutation_model, typename gamete_insertion_policy,
              typename gcont_t, typename mcont_t>
    std::size_t mutate_gamete_recycle(queue_type &recycling_bin,
                                      queue_type2 &gamete_recycling_bin,
                                      const gsl_rng *r, const double &mu,
                                      gcont_t &gametes, mcont_t &mutations,
                                      const size_t g,
                                      const mutation_model &mmodel,
                                      const gamete_insertion_policy &gpolicy);
}
#endif /* _FWDPP_MUTATION_HPP_ */
#include <fwdpp/mutation.tcc>
