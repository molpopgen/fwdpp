#ifndef FWDPP_GENETIC_MAP_UNIT_HPP
#define FWDPP_GENETIC_MAP_UNIT_HPP

#include <cstdint>
#include <functional>
#include <vector>
#include <gsl/gsl_rng.h>
#include <fwdpp/util/abstract_cloneable.hpp>

namespace fwdpp
{
    namespace detail
    {
        using genetic_map_unit_cast_function = std::function<double(double)>;

        inline genetic_map_unit_cast_function
        generate_genetic_map_unit_cast_function(bool discrete)
        {
            if (discrete)
                {
                    return [](double d) {
                        return static_cast<double>(static_cast<std::int64_t>(d));
                    };
                }

            return [](double d) { return d; };
        }
    }

    struct genetic_map_unit : public util::abstract_cloneable<genetic_map_unit>
    {
        genetic_map_unit() : util::abstract_cloneable<genetic_map_unit>()
        {
        }
        virtual void operator()(const gsl_rng *, std::vector<double> &) const = 0;
        virtual bool discrete() const = 0;
    };
} // namespace fwdpp

#endif
