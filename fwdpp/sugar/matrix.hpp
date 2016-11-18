#ifndef FWDPP_MATRIX_HPP_
#define FWDPP_MATRIX_HPP_

#include <utility>
#include <vector>
#include <type_traits>
#include <cstdint>
#include <algorithm>
#include <stdexcept>
#include <unordered_set>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/type_traits.hpp>
#include <fwdpp/sugar/poptypes/tags.hpp>

namespace KTfwd
{
    struct data_matrix
    {
        std::vector<std::int8_t> neutral, selected;
        std::vector<double> neutral_positions, selected_positions,
            neutral_popfreq, selected_popfreq;
        std::size_t nrow;
        data_matrix()
            : neutral{}, selected{}, neutral_positions{}, selected_positions{},
              neutral_popfreq{}, selected_popfreq{}, nrow{}
        {
        }
    };

    namespace data_matrix_details
    {
        template <typename gamete_t>
        void
        update_mutation_keys(std::unordered_set<std::size_t> &keys,
                             const typename gamete_t::mutation_container &a,
                             const std::vector<uint_t> &mcounts)
        {
            for (auto &&ai : a)
                {
                    if (mcounts[ai])
                        {
                            keys.insert(ai);
                        }
                }
        }

        template <typename dipvector_t, typename gcont_t>
        std::pair<std::vector<std::size_t>, std::vector<std::size_t>>
        mutation_keys(const dipvector_t &diploids,
                      const std::vector<std::size_t> &individuals,
                      const gcont_t &gametes,
                      const std::vector<uint_t> &mcounts,
                      const bool include_neutral, const bool include_selected,
                      const std::size_t, sugar::SINGLEPOP_TAG)
        {
            std::unordered_set<std::size_t> n, s;
            using gamete_t = typename gcont_t::value_type;
            for (auto &&ind : individuals)
                {
                    auto &dip = diploids[ind];
                    if (include_neutral)
                        {
                            update_mutation_keys<gamete_t>(
                                n, gametes[dip.first].mutations, mcounts);
                            update_mutation_keys<gamete_t>(
                                n, gametes[dip.second].mutations, mcounts);
                        }
                    if (include_selected)
                        {
                            update_mutation_keys<gamete_t>(
                                s, gametes[dip.first].smutations, mcounts);
                            update_mutation_keys<gamete_t>(
                                s, gametes[dip.second].smutations, mcounts);
                        }
                }
            return std::make_pair(
                std::vector<std::size_t>(n.begin(), n.end()),
                std::vector<std::size_t>(s.begin(), s.end()));
        }

        template <typename dipvector_t, typename gcont_t>
        std::pair<std::vector<std::size_t>, std::vector<std::size_t>>
        mutation_keys(const dipvector_t &diploids,
                      const std::vector<std::size_t> &individuals,
                      const gcont_t &gametes,
                      const std::vector<uint_t> &mcounts,
                      const bool include_neutral, const bool include_selected,
                      const std::size_t deme, sugar::METAPOP_TAG)
        {
            if (deme >= diploids.size())
                {
                    throw std::out_of_range("deme index out of range, "
                                            + std::string(__FILE__) + ' '
                                            + std::to_string(__LINE__));
                }
            // The i-th deme is a "singlepop", so we can use that function
            return mutation_keys(diploids[deme], individuals, gametes, mcounts,
                                 include_neutral, include_selected, deme,
                                 sugar::SINGLEPOP_TAG());
        }

        template <typename dipvector_t, typename gcont_t>
        std::pair<std::vector<std::size_t>, std::vector<std::size_t>>
        mutation_keys(const dipvector_t &diploids,
                      const std::vector<std::size_t> &individuals,
                      const gcont_t &gametes,
                      const std::vector<uint_t> &mcounts,
                      const bool include_neutral, const bool include_selected,
                      const std::size_t, sugar::MULTILOCPOP_TAG)
        {
            std::unordered_set<std::size_t> n, s;
            using gamete_t = typename gcont_t::value_type;
            for (auto &&ind : individuals)
                {
                    auto &dip = diploids[ind];
                    for (auto &&locus : dip)
                        {
                            if (include_neutral)
                                {
                                    update_mutation_keys<gamete_t>(
                                        n, gametes[locus.first].mutations,
                                        mcounts);
                                    update_mutation_keys<gamete_t>(
                                        n, gametes[locus.second].mutations,
                                        mcounts);
                                }
                            if (include_selected)
                                {
                                    update_mutation_keys<gamete_t>(
                                        s, gametes[locus.first].smutations,
                                        mcounts);
                                    update_mutation_keys<gamete_t>(
                                        s, gametes[locus.second].smutations,
                                        mcounts);
                                }
                        }
                }
            return std::make_pair(
                std::vector<std::size_t>(n.begin(), n.end()),
                std::vector<std::size_t>(s.begin(), s.end()));
        }

        template <typename poptype>
        data_matrix
        fill_data_matrix(const poptype &pop,
                         const std::vector<std::size_t> &individuals,
                         const std::vector<std::size_t> &neutral_keys,
                         const std::vector<std::size_t> &selected_keys,
                         const std::size_t deme)
        {
			data_matrix rv;

			return rv;
        }
    }

    template <typename poptype>
    std::pair<std::vector<std::size_t>, std::vector<std::size_t>>
    mutation_keys(const poptype &pop,
                  const std::vector<std::size_t> &individuals,
                  const bool include_neutral, const bool include_selected,
                  const std::size_t deme = 0)
    {
        return data_matrix_details::mutation_keys(
            pop.diploids, individuals, pop.gametes, pop.mcounts,
            include_neutral, include_selected, deme,
            typename poptype::popmodel_t());
    }

    template <typename poptype>
    data_matrix
    get_data_matrix(const poptype &pop,
                    const std::vector<std::size_t> &individuals,
                    const std::vector<std::size_t> &neutral_keys,
                    const std::vector<std::size_t> &selected_keys,
                    const std::size_t deme = 0)
    {
        return data_matrix_details::fill_data_matrix(
            pop, individuals, neutral_keys, selected_keys, deme);
    }
}

#endif
