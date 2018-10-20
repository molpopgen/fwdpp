#ifndef FWDPP_TS_MUTATE_TABLES_HPP
#define FWDPP_TS_MUTATE_TABLES_HPP


#include <gsl/gsl_randist.h>

#include <vector>
#include "definitions.hpp"
#include "mark_multiple_roots.hpp"
#include "table_collection.hpp"

namespace fwdpp
{
    namespace ts
    {
        template <typename rng, typename mfunction>
        unsigned
        mutate_tables(const rng &r, const mfunction &make_mutation,
                      table_collection &tables,
                      const std::vector<TS_NODE_INT> &samples, const double mu)
        {
            unsigned nmuts = 0;
            auto mr = mark_multiple_roots(tables, samples);
            for (auto &i : mr)
                {
                    auto dt = tables.node_table[i.first].time;
                    for (auto j : i.second)
                        {
                            double mean = dt * (j.second - j.first) * mu;
                            auto nm = gsl_ran_poisson(r.get(), mean);
                            nmuts += nm;
                            for (unsigned m = 0; m < nm; ++m)
                                {
                                    unsigned g = static_cast<unsigned>(
                                        gsl_ran_flat(r.get(), 1, dt + 1));
                                    auto k
                                        = make_mutation(j.first, j.second, g);
                                    tables.mutation_table.emplace_back(
                                        mutation_record{ i.first, k });
                                }
                        }
                }
            for (auto &e : tables.edge_table)
                {
                    auto ct = tables.node_table[e.child].time;
                    auto pt = tables.node_table[e.parent].time;
                    auto dt = ct - pt;
                    double mean = dt * (e.right - e.left) * mu;
                    auto nm = gsl_ran_poisson(r.get(), mean);
                    for (unsigned m = 0; m < nm; ++m)
                        {
                            unsigned g = static_cast<unsigned>(
                                gsl_ran_flat(r.get(), pt + 1, ct + 1));
                            auto k = make_mutation(e.left, e.right, g);
                            tables.mutation_table.emplace_back(
                                mutation_record{ e.child, k });
                        }
                    nmuts += nm;
                }
            return nmuts;
        }
    } // namespace ts
} // namespace fwdpp

#endif
