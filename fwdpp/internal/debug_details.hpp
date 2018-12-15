#ifndef FWDPP_INTERNAL_DEBUG_DETAILS_HPP
#define FWDPP_INTERNAL_DEBUG_DETAILS_HPP

#include <vector>
#include <algorithm>
#include <stdexcept>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/poptypes/tags.hpp>

/* We ignore unused variable warnings in this 
 * file.  In release mode, the functions are empty,
 * and will result in loads of unnecessary warnings
 */
#if defined(__clang__) && defined(__llvm__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"
#elif defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

namespace fwdpp
{
    namespace debug
    {
        namespace detail
        {
            template <typename gcont_t>
            void
            validate_sum_gamete_counts(const gcont_t &gametes,
                                       const uint_t expected_sum)
            {
#ifndef NDEBUG
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
#endif
            }

            template <typename mcont_t, typename iterator>
            void
            validate_mutation_key_ranges(const mcont_t &mutations,
                                         const iterator beg,
                                         const iterator end)
            {
#ifndef NDEBUG
                if (std::any_of(beg, end, [&mutations](const std::size_t m) {
                        return m >= mutations.size();
                    }))
                    {
                        throw std::runtime_error("FWDPP DEBUG: mutation key "
                                                 "out of range");
                    }
#endif
            }

            template <typename gamete_t>
            void
            gamete_is_extant(const gamete_t &gamete)
            {
#ifndef NDEBUG
                if (!gamete.n)
                    {
                        throw std::runtime_error(
                            "FWDPP DEBUG: unexpected extinct gamete");
                    }
#endif
            }

            template <typename gamete_t, typename mcont_t>
            void
            gamete_is_sorted(const gamete_t &g, const mcont_t &m)
            /*!
             * \brief Check that neutral mutation keys are sorted according to mutation
             * position
             */
            {
#ifndef NDEBUG
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
#endif
            }

            template <typename key_container, typename mcont_t>
            void
            gamete_data_valid(const key_container &keys,
                              const mcont_t &mutations,
                              const std::vector<uint_t> &mutcounts,
                              const bool expected_neutrality)
            {
#ifndef NDEBUG
                for (auto &k : keys)
                    {
                        if (!mutcounts[k])
                            {
                                throw std::runtime_error(
                                    "FWDPP DEBUG: extinct mutation in extant "
                                    "gamete");
                            }
                        if (mutations[k].neutral != expected_neutrality)
                            {
                                throw std::runtime_error(
                                    "FWDPP DEBUG: mutation neutrality field "
                                    "incorrect");
                            }
                    }
#endif
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
#ifndef NDEBUG
                detail::gamete_is_sorted(g, mutations);
                detail::gamete_data_valid(g.mutations, mutations, mutcounts,
                                          true);
                detail::gamete_data_valid(g.smutations, mutations, mutcounts,
                                          false);
#endif
            }

            template <typename diploid, typename gcont_t, typename mcont_t>
            void
            validate_pop_data_common(const diploid &dip,
                                     const gcont_t &gametes,
                                     const mcont_t &mutations,
                                     const std::vector<uint_t> &mutcounts)
            {
#ifndef NDEBUG
                if (!gametes[dip.first].n || !gametes[dip.second].n)
                    {
                        throw std::runtime_error(
                            "FWDPP DEBUG: gamete count is zero");
                    }
                gamete_is_sorted(gametes[dip.first], mutations);
                gamete_is_sorted(gametes[dip.second], mutations);
                gamete_data_valid(gametes[dip.first], mutations, mutcounts);
                gamete_data_valid(gametes[dip.second], mutations, mutcounts);
#endif
            }

            template <typename poptype>
            void
            validate_pop_data(const poptype &pop, poptypes::SINGLELOC_TAG)
            {
#ifndef NDEBUG
                for (const auto &d : pop.diploids)
                    {
                        validate_pop_data_common(d, pop.gametes, pop.mutations,
                                                 pop.mcounts);
                    }
#endif
            }

            template <typename poptype>
            void
            validate_pop_data(const poptype &pop, poptypes::MULTILOC_TAG)
            {
#ifndef NDEBUG
                for (const auto &d : pop.diploids)
                    {
                        for (const auto &locus : d)
                            {
                                validate_pop_data_common(locus, pop.gametes,
                                                         pop.mutations,
                                                         pop.mcounts);
                            }
                    }
#endif
            }

            template <typename poptype>
            void
            validate_pop_data(const poptype &pop)
            {
#ifndef NDEBUG
                validate_pop_data(pop, typename poptype::popmodel_t());
#endif
            }

            template <typename mutation_type>
            void
            check_mutation_neutrality(const mutation_type &mutation,
                                      const bool expected_neutrality)
            {
#ifndef NDEBUG
                if (mutation.neutral != expected_neutrality)
                    {
                        throw std::runtime_error(
                            "FWDPP DEBUG: mutation neutrality field "
                            "incorrect");
                    }
#endif
            }

            template <typename poptype>
            void
            all_gametes_extant(const poptype &pop,
                               const fwdpp::poptypes::SINGLELOC_TAG)
            {
#ifndef NEBUG
                for (auto &dip : pop.diploids)
                    {
                        if (pop.gametes[dip.first].n == 0
                            || pop.gametes[dip.second].n == 0)
                            {
                                throw std::runtime_error(
                                    "FWDPP DEBUG: diploid refers to "
                                    "extinct gamete");
                            }
                    }
#endif
            }

            template <typename poptype>
            void
            all_gametes_extant(const poptype &pop,
                               const fwdpp::poptypes::MULTILOC_TAG)
            {
#ifndef NEBUG
                for (auto &dip : pop.diploids)
                    {
                        for (auto &locus : dip)
                            {
                                if (pop.gametes[locus.first].n == 0
                                    || pop.gametes[locus.second].n == 0)
                                    {
                                        throw std::runtime_error(
                                            "FWDPP DEBUG: diploid refers to "
                                            "extinct gamete");
                                    }
                            }
                    }
#endif
            }
        } // namespace detail
    }     // namespace debug
} // namespace fwdpp

#if defined(__clang__) && defined(__llvm__)
#pragma clang diagnostic pop
#elif defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

#endif
