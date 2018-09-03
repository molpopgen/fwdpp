#ifndef FWDPP_SAMPLING_FUNCTIONS_DETAILS_HPP
#define FWDPP_SAMPLING_FUNCTIONS_DETAILS_HPP

#include <vector>
#include <cstdint>
#include <utility>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/data_matrix.hpp>

namespace fwdpp
{
    namespace fwdpp_internal
    {
        template <typename poptype>
        auto
        generate_filter_sort_keys(const poptype &pop,
                                  const std::vector<std::size_t> &individuals,
                                  const bool include_neutral,
                                  const bool include_selected,
                                  const bool remove_fixed)
            -> decltype(mutation_keys(pop, individuals, include_neutral,
                                      include_selected))
        {
            auto keys = mutation_keys(pop, individuals, include_neutral,
                                      include_selected);
            if (remove_fixed)
                {
                    const auto nsam = 2 * individuals.size();
                    const auto f
                        = [nsam](const std::pair<std::size_t, uint_t> &p) {
                              return p.second == nsam;
                          };
                    filter_keys(keys.first, f);
                    filter_keys(keys.second, f);
                }
            sort_keys(pop.mutations, keys.first);
            sort_keys(pop.mutations, keys.second);
            return keys;
        }
    }
}

#endif
