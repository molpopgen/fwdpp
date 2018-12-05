#ifndef FWDPP_SIM_PARAMS_HPP
#define FWDPP_SIM_PARAMS_HPP

#include <vector>
#include <utility>
#include <fwdpp/type_traits.hpp>
#include <fwdpp/forward_types.hpp>

namespace fwdpp
{
    template <typename genetic_value, typename mutation_function,
              typename recombination_function,
              typename interlocus_recombination_function>
    struct genetic_parameters
    {
        const genetic_value gvalue;
        const mutation_function generate_mutations;
        const recombination_function generate_breakpoints;
        const interlocus_recombination_function interlocus_recombination;
        std::vector<uint_t> neutral;
        std::vector<uint_t> selected;

        template <typename gv, typename mut, typename rec>
        genetic_parameters(gv&& gvalue_param, mut&& generate_mutations_param,
                           rec&& generate_breakpoints_param)
            : gvalue{ std::forward<gv>(gvalue_param) },
              generate_mutations{ std::forward<mut>(
                  generate_mutations_param) },
              generate_breakpoints{ std::forward<rec>(
                  generate_breakpoints_param) },
              interlocus_recombination{}, neutral{}, selected{}
        {
        }
    };

    template <typename genetic_value, typename mutation_function,
              typename recombination_function>
    inline genetic_parameters<genetic_value, mutation_function,
                              recombination_function, std::nullptr_t>
    make_genetic_parameters(
        genetic_value&& gvalue_param,
        mutation_function&& generate_mutations_param,
        recombination_function&& generate_breakpoints_param)
    {
        return genetic_parameters<genetic_value, mutation_function,
                                  recombination_function, std::nullptr_t>(
            std::move(gvalue_param), std::move(generate_mutations_param),
            std::move(generate_breakpoints_param));
    }

    template <typename genetic_value, typename mutation_function,
              typename recombination_function,
              typename interlocus_recombination>
    inline genetic_parameters<genetic_value, mutation_function,
                              recombination_function, interlocus_recombination>
    make_genetic_parameters(
        genetic_value&& gvalue_param,
        mutation_function&& generate_mutations_param,
        recombination_function&& generate_breakpoints_param,
        interlocus_recombination&& interlocus_rec)
    {
        return genetic_parameters<genetic_value, mutation_function,
                                  recombination_function, std::nullptr_t>(
            std::move(gvalue_param), std::move(generate_mutations_param),
            std::move(generate_breakpoints_param), std::move(interlocus_rec));
    }

} // namespace fwdpp

#endif

