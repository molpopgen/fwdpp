#ifndef FWDPP_TS_GENERATE_OFFSPRING_HPP
#define FWDPP_TS_GENERATE_OFFSPRING_HPP

#include <type_traits>
#include <gsl/gsl_rng.h>
#include <fwdpp/ts/get_parent_ids.hpp>
#include <fwdpp/ts/definitions.hpp>
#include <fwdpp/ts/table_collection.hpp>

namespace fwdpp
{
    namespace ts
    {
        struct ts_bookkeeper
        {
            double birth_time;
            std::int32_t offspring_deme;
            TS_NODE_INT first_parental_index, next_index;
        };

        template <typename breakpoint_function, typename new_mutation_fuction,
                  typename poptype, typename mrecbin, typename grecbin>
        inline TS_NODE_INT
        generate_offspring(const gsl_rng* r,
                           const std::pair<std::size_t, std::size_t> parents,
                           const ts_bookkeeper& bookkeeper,
                           const breakpoint_function& recombination_function,
                           const new_mutation_fuction& mutation_function,
                           poptype& pop,
                           typename poptype::diploid_t& offspring,
                           table_collection& tables,
                           mrecbin& mutation_recycling_bin,
                           grecbin& gamete_recycling_bin)
        {
            auto p1g1 = pop.diploids[parents.first].first;
            auto p1g2 = pop.diploids[parents.first].second;
            auto p2g1 = pop.diploids[parents.second].first;
            auto p2g2 = pop.diploids[parents.second].second;

            int swap1 = (gsl_rng_uniform(r) < 0.5) ? 1 : 0;
            int swap2 = (gsl_rng_uniform(r) < 0.5) ? 1 : 0;

            if (swap1)
                {
                    std::swap(p1g1, p1g2);
                }
            if (swap2)
                {
                    std::swap(p2g1, p2g2);
                }
            auto p1id = get_parent_ids(bookkeeper.first_parental_index,
                                       parents.first, swap1);
            auto p2id = get_parent_ids(bookkeeper.first_parental_index,
                                       parents.second, swap2);
            auto breakpoints = recombination_function();
        }
    } // namespace ts
} // namespace fwdpp

#endif
