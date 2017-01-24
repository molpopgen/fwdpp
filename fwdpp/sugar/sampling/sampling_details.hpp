#ifndef __FWDPP_SUGAR_SAMPLING_SAMPLING_DETAILS_HPP__
#define __FWDPP_SUGAR_SAMPLING_SAMPLING_DETAILS_HPP__
#include <vector>
#include <utility>
#include <string>
#include <algorithm>

namespace KTfwd
{

    namespace sugar
    {
        enum class treat_neutral
        {
            ALL,
            NEUTRAL,
            SELECTED
        };
    }

    template <typename vec_mutation_t>
    void
    add_fixations(sample_t &sample, const vec_mutation_t &fixations,
                  const unsigned nsam, const sugar::treat_neutral treat)
    {
        for (const auto &f : fixations)
            {
                if (treat == sugar::treat_neutral::ALL)
                    {
                        sample.emplace_back(f.pos, std::string(nsam, '1'));
                    }
                else if (treat == sugar::treat_neutral::NEUTRAL && f.neutral)
                    { // only add neutral mutations
                        sample.emplace_back(f.pos, std::string(nsam, '1'));
                    }
                else if (treat == sugar::treat_neutral::SELECTED && !f.neutral)
                    { // only add selected mutations
                        sample.emplace_back(f.pos, std::string(nsam, '1'));
                    }
            }
    }
    template <typename vec_mutation_t>
    void
    add_fixations(sample_t &sample, const vec_mutation_t &fixations,
                  const unsigned nsam, const sugar::treat_neutral treat,
                  const std::pair<double, double> &boundary)
    {
        for (const auto &f : fixations)
            {
                if (f.pos >= boundary.first && f.pos < boundary.second)
                    {
                        if (treat == sugar::treat_neutral::ALL)
                            {
                                sample.emplace_back(f.pos,
                                                    std::string(nsam, '1'));
                            }
                        else if (treat == sugar::treat_neutral::NEUTRAL
                                 && f.neutral)
                            { // only add neutral mutations
                                sample.emplace_back(f.pos,
                                                    std::string(nsam, '1'));
                            }
                        else if (treat == sugar::treat_neutral::SELECTED
                                 && !f.neutral)
                            { // only add selected mutations
                                sample.emplace_back(f.pos,
                                                    std::string(nsam, '1'));
                            }
                    }
            }
    }

    template <typename vec_mutation_t>
    void
    finish_sample(sample_t &sample, const vec_mutation_t &fixations,
                  const unsigned nsam, const bool removeFixed,
                  const sugar::treat_neutral treat)
    {
        if (!removeFixed)
            {
                add_fixations(sample, fixations, nsam, treat);
            }
        std::sort(sample.begin(), sample.end(),
                  [](const sample_site_t &a, const sample_site_t &b) noexcept {
                      return a.first < b.first;
                  });
    }
    template <typename vec_mutation_t>
    void
    finish_sample(sample_t &sample, const vec_mutation_t &fixations,
                  const unsigned nsam, const bool removeFixed,
                  const sugar::treat_neutral treat,
                  const std::pair<double, double> &boundary)
    {
        if (!removeFixed)
            {
                add_fixations(sample, fixations, nsam, treat, boundary);
            }
        std::sort(sample.begin(), sample.end(),
                  [](const sample_site_t &a, const sample_site_t &b) noexcept {
                      return a.first < b.first;
                  });
    }

    template <typename vec_mutation_t>
    void
    finish_sample(sep_sample_t &sample, const vec_mutation_t &fixations,
                  const unsigned nsam, const bool removeFixed,
                  const sugar::treat_neutral)
    {
        finish_sample(sample.first, fixations, nsam, removeFixed,
                      sugar::treat_neutral::NEUTRAL);
        finish_sample(sample.second, fixations, nsam, removeFixed,
                      sugar::treat_neutral::SELECTED);
    }

    // The next two are for multi-locus simulations
    template <typename vec_mutation_t>
    void
    finish_sample(
        std::vector<sample_t> &sample, const vec_mutation_t &fixations,
        const unsigned nsam, const bool removeFixed,
        const sugar::treat_neutral,
        const std::vector<std::pair<double, double>> &locus_boundaries)
    {
        for (std::size_t i = 0; i < sample.size(); ++i)
            {
                finish_sample(sample[i], fixations, nsam, removeFixed,
                              sugar::treat_neutral::ALL, locus_boundaries[i]);
            }
    }

    template <typename vec_mutation_t>
    void
    finish_sample(
        std::vector<sep_sample_t> &sample, const vec_mutation_t &fixations,
        const unsigned nsam, const bool removeFixed,
        const sugar::treat_neutral,
        const std::vector<std::pair<double, double>> &locus_boundaries)
    {
        for (std::size_t i = 0; i < sample.size(); ++i)
            {
                finish_sample(sample[i].first, fixations, nsam, removeFixed,
                              sugar::treat_neutral::NEUTRAL,
                              locus_boundaries[i]);
                finish_sample(sample[i].second, fixations, nsam, removeFixed,
                              sugar::treat_neutral::SELECTED,
                              locus_boundaries[i]);
            }
    }

    template <typename poptype, typename integer_type = std::size_t>
    sample_t
    sample_details(const poptype &p,
                   const std::vector<integer_type> &individuals,
                   const bool removeFixed)
    {
        sep_sample_t temp = fwdpp_internal::ms_sample_separate_single_deme(
            p.mutations, p.gametes, p.diploids, individuals,
            2 * individuals.size(), removeFixed);
        auto rv = std::move(temp.first);
        std::move(temp.second.begin(), temp.second.end(),
                  std::back_inserter(rv));
        finish_sample(rv, p.fixations, 2 * individuals.size(), removeFixed,
                      sugar::treat_neutral::ALL);
        return rv;
    }

    template <typename poptype, typename integer_type = std::size_t>
    std::vector<sample_t>
    sample_details(
        const poptype &p, const std::vector<integer_type> &individuals,
        const bool removeFixed,
        const std::vector<std::pair<double, double>> &locus_boundaries)
    {
        auto temp = fwdpp_internal::ms_sample_separate_mlocus(
            p.mutations, p.gamtes, p.diploids, individuals,
            2 * individuals.size(), removeFixed);
        if (!removeFixed && temp.size() != locus_boundaries.size())
            {
                throw std::runtime_error(
                    "incorrect number of elements in locus_boundaries");
            }
        std::vector<sample_t> rv;
        std::size_t j = 0;
        for (auto &i : temp)
            {
                rv.emplace_back(std::move(i.first));
                std::move(i.second.begin(), i.second.end(),
                          std::back_inserter(rv[j]));
            }
        finish_sample(rv, p.fixations, 2 * individuals.size(), removeFixed,
                      sugar::treat_neutral::ALL);
        return rv;
    }
}

#endif
