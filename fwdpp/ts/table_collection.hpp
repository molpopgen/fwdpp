#ifndef FWDPP_TS_TABLE_COLLECTION_HPP
#define FWDPP_TS_TABLE_COLLECTION_HPP

#include <cstdint>
#include "definitions.hpp"
#include "types/table_collection.hpp"

namespace fwdpp
{
    namespace ts
    {
        using table_collection = types::table_collection<std::int32_t>;
    }
}

#endif

