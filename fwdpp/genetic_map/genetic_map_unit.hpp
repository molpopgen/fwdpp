#ifndef FWDPP_GENETIC_MAP_UNIT_HPP
#define FWDPP_GENETIC_MAP_UNIT_HPP

#include <stdexcept>
#include <type_traits>
#include <vector>
#include <gsl/gsl_rng.h>
#include <fwdpp/util/abstract_cloneable.hpp>

/// \defgroup genetic_map_types Types to specify genetic maps
/// These types define a flexible API for building genetics maps
/// for simulations.
///
/// fwdpp::genetic_map_unit defines an interface class that is the
/// basis for a genetic map.
///
/// Multiple genetic map units may be collected in a fwdpp::genetic_map.
///
/// The library provides several subclasses of fwdpp::genetic_map_unit
/// that cover many common use cases.

namespace fwdpp
{
    namespace internal
    {
        template <typename T>
        inline void
        validate_minimum_interval_size(T b, T e, std::true_type)
        {
            if (e - b <= 1)
                {
                    throw std::invalid_argument(
                        "discrete intervals must have length greater than one.");
                }
        }

        template <typename T>
        inline void validate_minimum_interval_size(T, T, std::false_type)
        {
        }

        template <typename T>
        inline void
        validate_minimum_interval_size(T b, T e)
        {
            validate_minimum_interval_size(b, e, typename std::is_integral<T>::type());
        }
    }

    /// @struct genetic_map_unit genetic_map_unit.hpp fwdpp/genetic_map/genetic_map_unit.hpp
    /// @brief Interface class
    /// @ingroup genetic_map_types
    struct genetic_map_unit : public util::abstract_cloneable<genetic_map_unit>
    {
        /// Constructor
        genetic_map_unit() : util::abstract_cloneable<genetic_map_unit>()
        {
        }
        /// @note Future revisions may change the return type to void
        /// and allow for a reusable vector.
        virtual void operator()(const gsl_rng *, std::vector<double> &) const = 0;
    };
} // namespace fwdpp

#endif
