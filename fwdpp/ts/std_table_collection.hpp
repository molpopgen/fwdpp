#ifndef FWDPP_TS_STD_TABLE_COLLECTION_HPP
#define FWDPP_TS_STD_TABLE_COLLECTION_HPP

#warning("Deprecated: use <fwdpp/ts/table_collection.hpp> instead")

#include <vector>
#include "node.hpp"
#include "edge.hpp"
#include "site.hpp"
#include "mutation_record.hpp"
#include "detail/table_collection.hpp"

namespace fwdpp
{
    namespace ts
    {
        /// Alias for table collection backed by std::vector.
        /// \version 0.9.0

        using std_table_collection
            [[deprecated("use fwdpp::ts::table_collection from "
                         "<fwdpp/ts/table_collection.hpp> instead")]]
            = detail::table_collection<std::vector<node>, std::vector<edge>,
                                       std::vector<site>, std::vector<mutation_record>>;
    }
}

#endif
