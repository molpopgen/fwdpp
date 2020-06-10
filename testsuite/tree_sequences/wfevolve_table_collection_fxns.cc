#include "wfevolve_table_collection.hpp"

void
deaths_and_parents(const fwdpp::GSLrng_mt& rng, const std::vector<parent>& parents,
                   double psurvival, std::vector<birth>& births)
{
    births.clear();
    for (std::size_t i = 0; i < parents.size(); ++i)
        {
            if (gsl_rng_uniform(rng.get()) > psurvival)
                {
                    std::size_t parent0 = gsl_ran_flat(rng.get(), 0, parents.size());
                    std::size_t parent1 = gsl_ran_flat(rng.get(), 0, parents.size());
                    births.emplace_back(i, parents[parent0], parents[parent1]);
                }
        }
}

void
recombination_breakpoints(const fwdpp::GSLrng_mt& rng, double littler, double maxlen,
                          std::vector<double>& breakpoints)
{
    breakpoints.clear();
    auto nxovers = gsl_ran_poisson(rng.get(), littler);
    for (decltype(nxovers) i = 0; i < nxovers; ++i)
        {
            breakpoints.push_back(gsl_ran_flat(rng.get(), 0., maxlen));
        }
    std::sort(begin(breakpoints), end(breakpoints));
    if (!breakpoints.empty())
        {
            breakpoints.emplace_back(std::numeric_limits<double>::max());
        }
    return;
}
