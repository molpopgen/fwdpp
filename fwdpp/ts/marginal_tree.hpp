#pragma once

#include "types/marginal_tree.hpp"

namespace fwdpp
{
    namespace ts
    {
        template <typename SignedInteger>
        using marginal_tree = types::marginal_tree<SignedInteger>;

        template <typename SignedInteger>
        using sample_group_map = types::sample_group_map<SignedInteger>;
    }
}
