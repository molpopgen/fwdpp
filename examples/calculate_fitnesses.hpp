#ifndef FWDPP_EXAMPLES_CALCULATE_FITNESSES_HPP
#define FWDPP_EXAMPLES_CALCULATE_FITNESSES_HPP

#include <cstdint>
#include <vector>
#include <fwdpp/fitness_models.hpp>
#include <fwdpp/internal/gsl_discrete.hpp>

template <typename poptype, typename fitness_function>
inline fwdpp::fwdpp_internal::gsl_ran_discrete_t_ptr
calculate_fitnesses(poptype &pop, std::vector<double> &fitnesses,
                    const fitness_function &ff)
{
    auto N_curr = pop.diploids.size();
    fitnesses.resize(N_curr);
    for (std::size_t i = 0; i < N_curr; ++i)
        {
            fitnesses[i] = ff(pop.diploids[i], pop.gametes, pop.mutations);
        }
    auto lookup = fwdpp::fwdpp_internal::gsl_ran_discrete_t_ptr(
        gsl_ran_discrete_preproc(N_curr, &fitnesses[0]));
    return lookup;
}

#endif
