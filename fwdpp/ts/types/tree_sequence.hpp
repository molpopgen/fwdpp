#pragma once

#include <cstdint>
#include <memory>
#include <vector>
#include "table_collection.hpp"

namespace fwdpp
{
    namespace ts
    {
        namespace types
        {
            template <typename SignedInteger> class tree_sequence
            {
              private:
                std::shared_ptr<const types::table_collection<SignedInteger>> tables_;
                std::vector<SignedInteger> samples_;
                std::size_t num_trees_;

              public:
                tree_sequence(
                    std::shared_ptr<types::table_collection<SignedInteger>> tables)
                    : tables_{std::move(tables)}
                {
                }
            };
        }
    }
}
