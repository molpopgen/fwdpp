#ifndef FWDPP_TS_TABLE_COLLECTION_HPP
#define FWDPP_TS_TABLE_COLLECTION_HPP

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
        using table_collection
            = detail::table_collection<std::vector<node>, std::vector<edge>,
                                       std::vector<site>, std::vector<mutation_record>>;
    }
}

#endif

