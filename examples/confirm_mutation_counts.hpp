#ifndef FWDPP_EXAMPLES_CONFIRM_MUTATION_COUNTS
#define FWDPP_EXAMPLES_CONFIRM_MUTATION_COUNTS

#include <stdexcept>
#include <algorithm>
#include <vector>
#include <cstdint>
#include <fwdpp/internal/sample_diploid_helpers.hpp>
#include <fwdpp/ts/std_table_collection.hpp>

#ifndef NDEBUG
template <typename poptype>
void
confirm_mutation_counts(poptype &pop,
                        const fwdpp::ts::std_table_collection &tables)
{
    std::vector<std::size_t> keys;
    for (auto &mr : tables.mutations)
        {
            keys.push_back(mr.key);
        }
    std::sort(keys.begin(), keys.end());
    auto u = std::unique(keys.begin(), keys.end());
    if (u != keys.end())
        {
            throw std::runtime_error("redundant mutation keys");
        }

    for (auto &mr : tables.mutations)
        {
            if (mr.node == fwdpp::ts::TS_NULL_NODE)
                {
                    throw std::runtime_error(
                        "mutation node maps to null node");
                }

            if (tables.nodes[mr.node].time < pop.mutations[mr.key].g)
                {
                    throw std::runtime_error(
                        "mutation mapped to node pre-dating its origin");
                }
        }
    decltype(pop.mcounts) mc;
    fwdpp::fwdpp_internal::process_haploid_genomes(pop.haploid_genomes, pop.mutations, mc);
    for (std::size_t i = 0; i < pop.mcounts.size(); ++i)
        {
            if (!pop.mutations[i].neutral)
                {
                    assert(pop.mcounts[i] == mc[i]);
                }
        }
}
#else
template <typename poptype>
void
confirm_mutation_counts(poptype &, const fwdpp::ts::std_table_collection &)
{
}
#endif

#endif
