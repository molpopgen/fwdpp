#ifndef FWDPP_TS_VARIANT_FILLER_HPP
#define FWDPP_TS_VARIANT_FILLER_HPP

#include <cassert>
#include <vector>
#include <cstdint>
#include <algorithm>
#include "msprime_algo.hpp"

namespace fwdpp
{
    namespace ts
    {
        //TODO: rethink name.
        //variant_filler is a misnomer b/c what is really happening
        //is that samples are being tracked via top-2-bottom
        //traversal.
        //Best use of this type is probably composition
        //of a "visitor" for algorithmS.
        class variant_filler
        /// \brief Fill 0/1 genotype array for a mutation.
        {
          private:
            std::vector<TS_NODE_INT> sample_indexes;
            std::vector<TS_NODE_INT> node_stack;
            std::vector<std::int8_t> genotypes;
            TS_NODE_INT stack_top;

          public:
            variant_filler(const TS_NODE_INT nnodes,
                           const std::vector<TS_NODE_INT>& samples)
                : sample_indexes(nnodes, TS_NULL_NODE),
                  node_stack(nnodes + 1, TS_NULL_NODE),
                  genotypes(samples.size(), 0), stack_top(TS_NULL_NODE)
            {
                for (std::size_t i = 0; i < samples.size(); ++i)
                    {
                        auto s = samples[i];
                        if (s < 0 || s >= nnodes)
                            {
                                throw std::out_of_range("sample list contains "
                                                        "invalid node labels");
                            }
                        if (sample_indexes[s] != TS_NULL_NODE)
                            {
                                throw std::invalid_argument(
                                    "sample labels must be unique");
                            }
                        sample_indexes[static_cast<std::size_t>(s)]
                            = static_cast<TS_NODE_INT>(i);
                    }
            }

            inline std::uint32_t
            operator()(const marginal_tree& marginal, const TS_NODE_INT node)
            /// Returns true if node leads to any samples, false otherwise
            {
                //TODO need to compare this to an algorithmS-based version
                if (node < 0
                    || static_cast<std::size_t>(node) >= sample_indexes.size())
                    {
                        throw std::out_of_range("node value out of range");
                    }
                std::uint32_t variant_in_samples = 0;
                std::fill(std::begin(genotypes), std::end(genotypes), 0);
                const auto right = marginal.right_sample[node];
                for (auto i = marginal.left_sample[node]; i != TS_NULL_NODE;
                     i = marginal.next_sample[i])
                    {
                        auto sample_index = sample_indexes[i];
                        if (sample_index != TS_NULL_NODE)
                            {
                                ++variant_in_samples;
                                genotypes[sample_index] = 1;
                            }
                        if (i == right)
                            {
                                break;
                            }
                    }

                return variant_in_samples;
            }

            std::pair<const std::int8_t*, const std::int8_t*>
            view_genotypes() const
            {
                return std::make_pair(genotypes.data(),
                                      genotypes.data() + genotypes.size());
            }
        };
    } // namespace ts
} // namespace fwdpp

#endif
