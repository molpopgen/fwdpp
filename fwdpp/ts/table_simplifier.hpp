#ifndef FWDPP_TS_TABLE_SIMPLIFIER_HPP
#define FWDPP_TS_TABLE_SIMPLIFIER_HPP

#include <cmath>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cstddef>
#include <tuple>
#include <stdexcept>
#include "definitions.hpp"
#include "simplification/simplification.hpp"
#include "simplify_tables.hpp"

namespace fwdpp
{
    namespace ts
    {
        template <typename TableCollectionType> class table_simplifier
        /*! \brief Implements the simplification algorithm of \cite Kelleher2018-fu
         *
         *  Keeping a persistent simplifier around during a simplification avoids many 
         *  repeated large memory allocations.  On a single-core machine, the more often
         *  you simplify, the more this matters.  When running many simulations on a many-core 
         *  machine, keeping a persistent object reduces competition between threads for big
         *  memory chunks, which should also lead to a performance boost.
         *
         *  Many of the implementation details are private functions, which are subject to change
         *  without notice.
         *
         *  \version 0.7.0 Added to fwdpp.
         *  \version 0.8.0 Added ancestry_list, which greatly reduces memory overhead.
         *  \version 0.8.0 Remove need for temporary output edge table.
         *  \version 0.9.0 Added typename TableCollectionType and refactored as a wrapper around
         *                 standalone functions.
         */
        {
          private:
            using simplifier_internal_state
                = simplification::simplifier_internal_state<TableCollectionType>;

            simplifier_internal_state _state;

          public:
            table_simplifier() : _state{}
            {
            }

            std::pair<std::vector<table_index_t>, std::vector<std::size_t>>
            simplify(TableCollectionType& tables,
                     const std::vector<table_index_t>& samples)
            /// Simplify algorithm is approximately the same
            /// logic as used in msprime 0.6.0
            ///
            /// \param tables A table_collection
            /// \param samples A list of sample (node) ids.
            /// \version 0.7.1 Throw exception if a sample is recorded twice
            /// \version 0.7.3 Return value is now a pair containing the
            /// node ID map and a vector of keys to mutations preserved in
            /// mutation tables
            {
                simplify_tables_output simplification_output;
                simplify_tables(samples, simplification_flags{}, _state, tables,
                                simplification_output);

                return std::make_pair(
                    std::move(simplification_output.idmap),
                    std::move(simplification_output.preserved_mutations));
            }
        };

        template <typename TableCollectionType>
        inline table_simplifier<TableCollectionType>
        make_table_simplifier(const TableCollectionType&)
        /// Convenience function to generate a simplifier.
        /// \version 0.9.0 Added to library
        {
            return table_simplifier<TableCollectionType>();
        }
    } // namespace ts
} // namespace fwdpp

#endif
