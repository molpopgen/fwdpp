#ifndef FWDPP_TS_GENERATE_DATA_MATRIX_HPP
#define FWDPP_TS_GENERATE_DATA_MATRIX_HPP

#include <vector>
#include <type_traits>
#include <cstdint>
#include <stdexcept>
#include <fwdpp/data_matrix.hpp>
#include "table_collection.hpp"
#include "visit_sites.hpp"
#include "exceptions.hpp"
#include "marginal_tree_functions/samples.hpp"
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
        /// \version 0.8.0 No longer requires mutation vector. Function body re-implemented.
        {
            if (!(stop > start))
                {
                    throw std::invalid_argument("invalid position range");
                }
            std::vector<std::int8_t> genotypes(samples.size());
            detail::data_matrix_filler visitor(samples.size(), record_neutral,
                                               record_selected, skip_fixed);
            visit_sites(tables, std::forward<SAMPLES>(samples),
                        std::ref(visitor), start, stop);

            return visitor.data();
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
