#ifndef FWDPP_TS_MAKE_SIMPLIFIER_STATE_HPP
#define FWDPP_TS_MAKE_SIMPLIFIER_STATE_HPP

#include <fwdpp/ts/simplification/simplification.hpp>

namespace fwdpp
{
    namespace ts
    {
        template <typename TableCollectionType>
        inline simplification::simplifier_internal_state<TableCollectionType>
        make_simplifier_state(const TableCollectionType& tables)
        {
            return simplification::make_simplifier_internal_state(tables);
        }
    }
}

#endif
