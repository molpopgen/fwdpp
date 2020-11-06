#ifndef FWDPP_GENETIC_MAP_UNIT_HPP
#define FWDPP_GENETIC_MAP_UNIT_HPP

#include <vector>
#include <gsl/gsl_rng.h>
#include <fwdpp/util/abstract_cloneable.hpp>

/// \defgroup genetic_map_unit_types Types to specify genetic maps
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
    struct genetic_map_unit : public util::abstract_cloneable<genetic_map_unit>
    /// \brief Interface class
    /// \ingroup genetic_map_unit_types
    {
        genetic_map_unit() : util::abstract_cloneable<genetic_map_unit>()
        {
        }
        /// \note Future revisions may change the return type to void
        /// and allow for a reusable vector.
        virtual void operator()(const gsl_rng *, std::vector<double> &) const = 0;
    };
} // namespace fwdpp

#endif
