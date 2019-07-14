#ifndef FWDPP_TS_GENERATE_DATA_MATRIX_HPP
#define FWDPP_TS_GENERATE_DATA_MATRIX_HPP

#include <vector>
#include <type_traits>
#include <cstdint>
#include <stdexcept>
#include <fwdpp/data_matrix.hpp>
#include "table_collection.hpp"
#include "tree_visitor.hpp"
#include "detail/generate_data_matrix_details.hpp"

namespace fwdpp
{
    namespace ts
    {
        template <typename SAMPLES>
        data_matrix
        generate_data_matrix(const table_collection& tables, SAMPLES&& samples,
                             const bool record_neutral,
                             const bool record_selected, const bool skip_fixed,
                             const double start, const double stop)
        /// \todo Document
        /// \version 0.7.0 Added to library
        /// \version 0.7.1 Change behavior to skip sites fixed in the sample
        /// \version 0.7.4 Add [start, stop) arguments. Add option to skip fixed variants.
        /// \version 0.8.0 No longer requires mutation vector.
        {
            if (!(stop > start))
                {
                    throw std::invalid_argument("invalid position range");
                }
            auto mut = tables.mutation_table.cbegin();
            const auto mut_end = tables.mutation_table.cend();
            auto site_table_location = begin(tables.site_table);
            const auto site_table_end = end(tables.site_table);
            tree_visitor tv(tables, std::forward<SAMPLES>(samples),
                            update_samples_list(true));
            detail::data_matrix_filler processor(
                tables.site_table.begin(), record_neutral, record_selected,
                skip_fixed, samples.size());
            while (tv())
                {
                    const auto& tree = tv.tree();
                    // Advance site table to start of tree:
                    while (site_table_location < site_table_end
                           && site_table_location->position < tree.left)
                        {
                            ++site_table_location;
                        }
                    // Process mutations on this tree
                    for (; site_table_location < site_table_end
                           && site_table_location->position < tree.right;
                         ++site_table_location)
                        {
                            auto pos = site_table_location->position;
                            if (pos >= start)
                                {
                                    if (!(pos < stop))
                                        {
                                            return processor.data();
                                        }
                                    //catch the mutation table up
                                    //to the site table
                                    while (mut < mut_end
                                           && tables.site_table[mut->site]
                                                      .position
                                                  != pos)
                                        {
                                            ++mut;
                                        }
                                    mut = processor(tree, *site_table_location,
                                                    mut, mut_end);
                                }
                        }
                    // Terminate early if we're lucky
                    if (!(mut < mut_end))
                        {
                            break;
                        }
                }
            return processor.data();
        }

        template <typename SAMPLES>
        data_matrix
        generate_data_matrix(const table_collection& tables, SAMPLES&& samples,
                             const bool record_neutral,
                             const bool record_selected, const bool skip_fixed)
        {
            return generate_data_matrix(
                tables, std::forward<SAMPLES>(samples), record_neutral,
                record_selected, skip_fixed, 0., tables.genome_length());
        }

    } // namespace ts
} // namespace fwdpp

#endif
