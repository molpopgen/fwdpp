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
    /// @class genetic_map genetic.hpp fwdpp/genetic_map/genetic_map.hpp
    /// @brief A genetic map is a container of fwdpp::genetic_map_unit objects.
    /// @ingroup genetic_map_types
    /// @headerfile genetic_map.hpp <fwdpp/genetic_map/genetic_map.hpp>
    class genetic_map
    {
      private:
        std::vector<std::unique_ptr<genetic_map_unit>> callbacks;

      public:
        genetic_map(std::vector<std::unique_ptr<genetic_map_unit>>&& c)
            : callbacks(std::move(c))
        /// Constructor
        /// @param c ``rvalue`` reference to vector of genetic_map_unit objects.
        {
        }

        genetic_map() : callbacks{}
        ///
        /// Construct an empty genetic_map
        {
        }

        void
        add_callback(std::unique_ptr<genetic_map_unit>&& gu)
        /// Add a new callback by moving the input
        {
            callbacks.emplace_back(std::move(gu));
        }

        void
        add_callback(const genetic_map_unit& gu)
        /// Add a new callback by cloning the input
        {
            callbacks.emplace_back(gu.clone());
        }

        std::vector<double>
        operator()(const gsl_rng* r) const
        /// @note Future revisions may change the return type to void
        /// and allow for a reusable vector.
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
        /// Return number of stored callbacks
        {
            return callbacks.size();
        }
    };
} // namespace fwdpp

#endif
