#ifndef FWDPP_GENETIC_MAP_HPP__
#define FWDPP_GENETIC_MAP_HPP__

#include <memory>
#include <vector>
#include <utility>
#include <limits>
#include <algorithm>
#include "genetic_map_unit.hpp"

namespace fwdpp
{
    class genetic_map
    {
      private:
        std::vector<std::unique_ptr<genetic_map_unit>> callbacks;

      public:
        genetic_map(std::vector<std::unique_ptr<genetic_map_unit>>&& c)
            : callbacks(std::move(c))
        {
        }
        genetic_map() : callbacks{} {}

        void
        add_callback(const genetic_map_unit& gu)
        {
            callbacks.emplace_back(gu.clone());
        }

        std::vector<double>
        operator()(const gsl_rng* r) const
        {
            std::vector<double> breakpoints;
            for (auto& c : callbacks)
                {
                    c->operator()(r, breakpoints);
                }
            std::sort(begin(breakpoints), end(breakpoints));
            if (!breakpoints.empty())
                {
                    breakpoints.push_back(std::numeric_limits<double>::max());
                }
            return breakpoints;
        }

        std::size_t
        size() const
        {
            return callbacks.size();
        }
    };
} // namespace fwdpp

#endif
