#ifndef FWDPP_TS_MAKE_SIMPLIFIER_STATE_HPP
#define FWDPP_TS_MAKE_SIMPLIFIER_STATE_HPP

#include <fwdpp/ts/simplification/simplification.hpp>
#include "types/table_collection.hpp"

namespace fwdpp
{
    namespace ts
    {
        template <typename SignedInteger>
        inline simplification::simplifier_internal_state<SignedInteger>
        make_simplifier_state(const types::table_collection<SignedInteger>& tables)
        {
            return simplification::make_simplifier_internal_state(tables);
        }
    }
}

#endif
