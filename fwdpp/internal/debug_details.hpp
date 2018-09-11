#ifndef FWDPP_INTERNAL_DEBUG_DETAILS_HPP
#define FWDPP_INTERNAL_DEBUG_DETAILS_HPP

#include <vector>
#include <algorithm>
#include <stdexcept>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/sugar/poptypes/tags.hpp>

namespace fwdpp
{
    namespace debug
    {
        namespace detail
        {
#ifndef NDEBUG
            template <typename gcont_t>
            void
            validate_sum_gamete_counts(const gcont_t &gametes,
                                       const uint_t expected_sum)
            {
                uint_t s = 0;
                for (auto &g : gametes)
                    {
                        s += g.n;
                    }
                if (s != expected_sum)
                    {
                        throw std::runtime_error(
                            "FWDPP DEBUG: unexpectd sum of gamete counts");
                    }
            }

            template <typename mcont_t, typename iterator>
            void
            validate_mutation_key_ranges(const mcont_t &mutations,
                                         const iterator beg,
                                         const iterator end)
            {
                if (std::any_of(beg, end, [&mutations](const std::size_t m) {
                        return m >= mutations.size();
                    }))
                    {
                        throw std::runtime_error("FWDPP DEBUG: mutation key "
                                                 "out of range");
                    }
            }

            template <typename gamete_t>
            void
            gamete_is_extant(const gamete_t &gamete)
            {
                if (!gamete.n)
                    {
                        throw std::runtime_error("FWDPP DEBUG: unexpected extinct gamete");
                    }
            }

            template <typename gamete_t, typename mcont_t>
            void
            gamete_is_sorted(const gamete_t &g, const mcont_t &m)
            /*!
             * \brief Check that neutral mutation keys are sorted according to mutation
             * position
             */
            {
                const auto comp = [&m](const size_t i, const size_t j) {
                    return m[i].pos <= m[j].pos;
                };

                if (!std::is_sorted(g.mutations.begin(), g.mutations.end(),
                                    comp))
                    {
                        throw std::runtime_error(
                            "FWDPP DEBUG: neutral mutation keys not sorted");
                    }
                if (!std::is_sorted(g.smutations.begin(), g.smutations.end(),
                                    comp))
                    {
                        throw std::runtime_error(
                            "FWDPP DEBUG: selected mutation keys not sorted");
                    }
            }

            template <typename key_container, typename mcont_t>
            void
            gamete_data_valid(const key_container &keys,
                              const mcont_t &mutations,
                              const std::vector<uint_t> &mutcounts,
                              const bool expected_neutrality)
            {
                for (auto &k : keys)
                    {
                        if (!mutcounts[k])
                            {
                                throw std::runtime_error(
                                    "FWDPP DEBUG: extinct mutation in extant gamete");
                            }
                        if (mutations[k].neutral != expected_neutrality)
                            {
                                throw std::runtime_error(
                                    "FWDPP DEBUG: mutation neutrality field incorrect");
                            }
                    }
            }

            template <typename gamete_t, typename mcont_t>
            void
            gamete_data_valid(const gamete_t &g, const mcont_t &mutations,
                              const std::vector<uint_t> &mutcounts)
            /*
      \brief Check that "neutral" and "non-neutral" mutations are where we
      expect them to be.
     */
            {
                detail::gamete_is_sorted(g, mutations);
                detail::gamete_data_valid(g.mutations, mutations, mutcounts,
                                          true);
                detail::gamete_data_valid(g.smutations, mutations, mutcounts,
                                          false);
            }

            template <typename diploid, typename gcont_t, typename mcont_t>
            void
            validate_pop_data_common(const diploid &dip,
                                     const gcont_t &gametes,
                                     const mcont_t &mutations,
                                     const std::vector<uint_t> &mutcounts)
            {
                if (!gametes[dip.first].n || !gametes[dip.second].n)
                    {
                        throw std::runtime_error("FWDPP DEBUG: gamete count is zero");
                    }
                gamete_is_sorted(gametes[dip.first], mutations);
                gamete_is_sorted(gametes[dip.second], mutations);
                gamete_data_valid(gametes[dip.first], mutations, mutcounts);
                gamete_data_valid(gametes[dip.second], mutations, mutcounts);
            }

            template <typename poptype>
            void
            validate_pop_data(const poptype &pop, sugar::SINGLELOC_TAG)
            {
                for (const auto &d : pop.diploids)
                    {
                        validate_pop_data_common(d, pop.gametes, pop.mutations,
                                                 pop.mcounts);
                    }
            }

            template <typename poptype>
            void
            validate_pop_data(const poptype &pop, sugar::MULTILOC_TAG)
            {
                for (const auto &d : pop.diploids)
                    {
                        for (const auto &locus : d)
                            {
                                validate_pop_data_common(locus, pop.gametes,
                                                         pop.mutations,
                                                         pop.mcounts);
                            }
                    }
            }

            template <typename poptype>
            void
            validate_pop_data(const poptype &pop)
            {
                validate_pop_data(pop, typename poptype::popmodel_t());
            }
#else
            template <typename gcont_t>
            void
            validate_sum_gamete_counts(const gcont_t & /*gametes*/,
                                       const uint_t /*expected_sum*/)
            {
            }

            template <typename mcont_t, typename iterator>
            void
            validate_mutation_key_ranges(const mcont_t & /*mutations*/,
                                         const iterator /*beg*/,
                                         const iterator /*end*/)
            {
            }

            template <typename gamete_t>
            void
            gamete_is_extant(const gamete_t & /*gamete*/)
            {
            }

            template <typename gamete_t, typename mcont_t>
            void
            gamete_is_sorted(const gamete_t &, const mcont_t &)
            /*!
              \brief Check that neutral mutation keys are sorted according to mutation
              position
            */
            {
            }

            template <typename gamete_t, typename mcont_t>
            void
            gamete_data_valid(const gamete_t & /*gamete*/,
                              const mcont_t & /*mutations*/,
                              const std::vector<uint_t> & /*mutcounts*/)
            {
            }

            template <typename poptype>
            void
            validate_pop_data(const poptype & /*pop*/)
            {
            }
#endif
        } // namespace detail
    }     // namespace debug
} // namespace fwdpp

#endif
