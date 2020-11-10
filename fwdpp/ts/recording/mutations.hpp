#ifndef FWDPP_TS_RECORDING_MUTATIONS_HPP
#define FWDPP_TS_RECORDING_MUTATIONS_HPP

#include <limits>
#include <vector>
#include <stdexcept>
#include <fwdpp/ts/definitions.hpp>

namespace fwdpp
{
    namespace ts
    {
        template <typename TableCollectionType, typename MutationContainerType>
        void
        record_mutations_infinite_sites(
            const table_index_t u, const MutationContainerType& mutations,
            const std::vector<std::uint32_t>& new_mutation_keys,
            TableCollectionType& tables)
        /// \version Added in 0.8.0
        /// FIXME: move this to a tree seq recording header
        {
            std::int8_t ancestral_state = 0, derived_state = 1;
            for (auto& k : new_mutation_keys)
                {
                    auto site
                        = tables.emplace_back_site(mutations[k].pos, ancestral_state);
                    if (site >= std::numeric_limits<table_index_t>::max())
                        {
                            throw std::invalid_argument("site index out of range");
                        }
                    tables.push_back_mutation(u, k, site, derived_state,
                                              mutations[k].neutral);
                }
        }
    }
}

#endif
