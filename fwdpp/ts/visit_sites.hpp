#ifndef FWDPP_TS_VISIT_SITES_HPP
#define FWDPP_TS_VISIT_SITES_HPP

#include "tree_visitor.hpp"

namespace fwdpp
{
    namespace ts
    {
        template <typename TableCollectionType, typename F, typename SAMPLES>
        void
        visit_sites(const TableCollectionType& tables, SAMPLES&& samples, const F& f,
                    const double from, const double to)
        /// \brief Apply a function to all sites in range [from, to)
        /// \param tables A ts::table_collection
        /// \param samples A list of samples
        /// \param f A function.  See below
        /// \param from Start of range (inclusive)
        /// \param to End of range (exclusive)
        ///
        /// The type of \a f must be equavalent to:
        ///
        /// \code{.cpp}
        /// std::function<void(const ts::marginal_tree&, const ts::site&,
        ///     typename TableCollectionType::mutation_table::const_iterator
        ///     typename TableCollectionType::mutation_table::const_iterator)>
        /// \endcode
        ///
        /// The last two iterators define the range of mutations in \tables correpsonding
        /// to the ts::site object passed in.  The range is defined following standard
        /// C++ conventions for bidirectional iterators from a vector.
        /// \version 0.8.0 Added to library
        {
            tree_visitor<TableCollectionType> tv(tables, std::forward<SAMPLES>(samples),
                                                 update_samples_list(true));
            auto current_site = begin(tables.sites);
            auto current_mutation = begin(tables.mutations);
            using site = typename TableCollectionType::site_table::value_type;
            using mutation_record =
                typename TableCollectionType::mutation_table::value_type;
            while (tv())
                {
                    if (tv.tree().left >= from && tv.tree().left < to)
                        {
                            if (current_site->position < tv.tree().left)
                                {
                                    current_site = std::lower_bound(
                                        current_site, end(tables.sites), tv.tree().left,
                                        [](const site& s, const double p) {
                                            return s.position < p;
                                        });
                                    if (current_site < end(tables.sites))
                                        {
                                            current_mutation = std::lower_bound(
                                                begin(tables.mutations),
                                                end(tables.mutations),
                                                current_site->position,
                                                [&tables](const mutation_record& m,
                                                          const double p) {
                                                    return tables.sites[m.site].position
                                                           < p;
                                                });
                                        }
                                }
                            while (current_site < end(tables.sites)
                                   && current_site->position < tv.tree().right
                                   && current_site->position < to)
                                {
                                    while (
                                        current_mutation < end(tables.mutations)
                                        && tables.sites[current_mutation->site].position
                                               < current_site->position)
                                        {
                                            ++current_mutation;
                                        }
                                    auto m = std::find_if(
                                        current_mutation, end(tables.mutations),
                                        [&tables,
                                         current_site](const mutation_record& mr) {
                                            return tables.sites[mr.site].position
                                                   != current_site->position;
                                        });
                                    f(tv.tree(), *current_site, current_mutation, m);
                                    current_mutation = m;
                                    ++current_site;
                                }
                        }
                    if (current_site == end(tables.sites) || current_site->position >= to
                        || tv.tree().left >= to)
                        {
                            return;
                        }
                }
        }
    } // namespace ts
} // namespace fwdpp

#endif
