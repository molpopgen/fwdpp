#ifndef FWDPP_MUTATE_RECOMBINE_HPP__
#define FWDPP_MUTATE_RECOMBINE_HPP__

#include <vector>
#include <algorithm>
#include <fwdpp/forward_types.hpp>
#include <fdwpp/internal/mutation_internal.hpp>

namespace KTfwd
{
    template <typename queue_type, typename queue_type2,
              typename mutation_model, typename gcont_t, typename mcont_t>
    std::vector<uint_t>
    generate_new_mutations(queue_type &recycling_bin, const gsl_rng *r,
                           const double &mu, 
                           mcont_t &mutations, const std::size_t g,
                           const mutation_model &mmodel)
    /// Return a vector of keys to new mutations.  The keys
    /// will be sorted according to mutation postition.
    {
        unsigned nm = gsl_ran_poisson(r, mu);
        std::vector<uint_t> rv;
        rv.reserve(nm);
        for (unsigned i = 0; i < nm; ++i)
            {
                rv.emplace_back(fwdpp_internal::mmodel_dispatcher(
                    mmodel, g, mutations, recycling_bin));
            }
        std::sort(mutations.begin(), mutations.end(),
                  [&mutations](const uint_t a, const uint_t b) {
                      return mutations[a].pos < mutations[b].pos;
                  });
        return rv;
    }
}

#endif
