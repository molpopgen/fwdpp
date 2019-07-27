#ifndef FWDPP_GENETIC_MAP_UNIT_HPP
#define FWDPP_GENETIC_MAP_UNIT_HPP

#include <vector>
#include <gsl/gsl_rng.h>
#include <fwdpp/util/abstract_cloneable.hpp>

namespace fwdpp
{
    struct genetic_map_unit : public util::abstract_cloneable<genetic_map_unit>
    {
        genetic_map_unit() : util::abstract_cloneable<genetic_map_unit>() {}
        virtual void operator()(const gsl_rng *,
                                std::vector<double> &) const = 0;
    };
} // namespace fwdpp

#endif
