#ifndef FWDPP_TS_MUTATE_TABLES_HPP
#define FWDPP_TS_MUTATE_TABLES_HPP

#include <gsl/gsl_randist.h>

#include <vector>
#include "definitions.hpp"
#include "mark_multiple_roots.hpp"
#include "mutation_tools.hpp"
#include "table_collection_functions.hpp"

namespace fwdpp
{
    namespace ts
    {
        template <typename TableCollectionType, typename Samples, typename rng,
                  typename mfunction>
        unsigned
        mutate_tables(const rng &r, const mfunction &make_mutation,
                      TableCollectionType &tables, Samples &&samples, const double mu)
        /// \brief Apply a mutation scheme to add neutral mutations to a fwdpp::ts::table_collection.
        ///
        /// \param r fwdpp::GSLrng_t
        /// \param make_mutation A mutation function.  See below.
        /// \param tables A fwdpp::ts::table_collection
        /// \param samples A list of sample nodes corresponding to "currently-alive" nodes
        /// \param mu Mutation rate (per gamete, per generation)
        ///
        /// \version 0.7.0 Added to library
        /// \version 0.7.2 Return immediately if mutation rate is not > 0
        /// \version 0.9.0 Added typename TableCollectionType
        ///
        /// This function uses the edge_table to apply mutations to all parent/child connections.
        /// The mutations should be neutral, although that requirement is not enforced. (It simply
        /// makes little sense to apply selected mutations to the entire edge table post-hoc.)
        ///
        /// The result of this function is to populate \a tables.mutations with neutral variants.
        ///
        /// The parameter \a make_mutation is a function that must conform to
        /// std::function<new_variant_record(double, double, double)>.  The three arguments
        /// are interpreted as "left", "right", and "time".  The function's return value
        /// allows us to properly update the site and mutation tables in \a tables.
        /// The new mutation must have a position on the half-open interval [left, right). The
        /// "time" parameger represents the generation when the mutation arose and will be
        /// uniformly assigned between the parental birth time and that of the child.
        /// This mutation's origin time will be on the half-open interval (parental birth time,
        /// child birth time].
        ///
        /// The preceding paragraph implies that \a make_mutation is repsonsible for required operations
        /// related to mutation recycling, etc..  The lambda "neutral_variant_maker" in the example
        /// program wfts.cc will probably be of help here.
        ///
        /// It is common that, at the end of a forward-time simulation, that not all trees are
        /// completely coalesced.  This happens because the distribution on the time to MRCA
        /// has very high variance and, due to recombination, some fraction of trees will not be
        /// coalesced after, say, 10N generations of evolving a single Wright-Fisher deme.  Further,
        /// the simlification algorithm will "push" the most ancient nodes on these "marginal forests"
        /// forwards in time to the most recent ancestral node of each tree in the forest.  Thus, naively
        /// mutating the edge table will place too few mutations on these parts of the genome.  This
        /// function corrects for the presence of marginal forests, as described below.
        ///
        /// The \a samples list is used to identify "marginal forests", meaning marginal trees
        /// that are not completely coalesced. (The finding is done via a call to
        /// fwdpp::ts::mark_multiple_roots.)  The MRCA nodes on trees in marginal forests have additional
        /// mutations on them, representing evolution from time point zero (the beginning of the simulation)
        /// to the MRCA node's time.
        ///
        /// Note that two alternatives exist that will render the treatment of marginal forests unnecessary:
        /// 1. Simulate for longer.  For example, CDF of TMRCA under Wright-Fisher is quite close to 1 at ~20N generations.
        /// Thus, simulating longer means that fewer and fewer marginal trees will be forests.
        /// 2. Start the simulation with an existing, completely-coalesced, tree sequence. As far as fwdpp is concerned,
        /// this is left as an "exercise for the reader" at the moment.
        {
            unsigned nmuts = 0;
            if (!(mu > 0.0))
                {
                    return nmuts;
                }
            auto mr = mark_multiple_roots(tables, std::forward<Samples>(samples));
            const double L = tables.genome_length();
            const auto make_origin_time
                = [&r](const double parent_time, const double child_time) {
                      auto fc = std::floor(child_time);
                      auto cp = std::ceil(parent_time);
                      auto dt = fc - cp;
                      if (dt <= 0.0)
                          {
                              if (fc <= parent_time)
                                  {
                                      return std::numeric_limits<double>::quiet_NaN();
                                  }
                              return fc;
                          }
                      auto x = gsl_ran_flat(r.get(), 0, dt);
                      if (fc - x <= parent_time)
                          {
                              return std::numeric_limits<double>::quiet_NaN();
                          }
                      return fc - x;
                  };
            for (auto &i : mr)
                {
                    auto dt = tables.nodes[i.first].time;
                    for (auto j : i.second)
                        {
                            double mean = dt * (j.second - j.first) * mu / L;
                            auto nm = gsl_ran_poisson(r.get(), mean);
                            nmuts += nm;
                            for (unsigned m = 0; m < nm; ++m)
                                {
                                    auto g = make_origin_time(0, dt);
                                    new_variant_record r
                                        = make_mutation(j.first, j.second, g);
                                    auto newsite = tables.emplace_back_site(
                                        r.s.position, r.s.ancestral_state);
                                    tables.emplace_back_mutation(i.first, r.key, newsite,
                                                                 r.derived_state,
                                                                 r.neutral);
                                }
                        }
                }
            for (auto &e : tables.edges)
                {
                    auto ct = tables.nodes[e.child].time;
                    auto pt = tables.nodes[e.parent].time;
                    auto dt = ct - pt;
                    double mean = dt * (e.right - e.left) * mu / L;
                    auto nm = gsl_ran_poisson(r.get(), mean);
                    for (unsigned m = 0; m < nm; ++m)
                        {
                            auto g = make_origin_time(pt, ct);
                            new_variant_record r = make_mutation(e.left, e.right, g);
                            auto site = tables.emplace_back_site(r.s.position,
                                                                 r.s.ancestral_state);
                            tables.emplace_back_mutation(e.child, r.key, site,
                                                         r.derived_state, r.neutral);
                        }
                    nmuts += nm;
                }
            sort_mutation_table_and_rebuild_site_table(tables);
            return nmuts;
        }
    } // namespace ts
} // namespace fwdpp

#endif
