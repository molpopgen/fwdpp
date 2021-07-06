#ifndef TS_EXAMPLES_TYPES_HPP
#define TS_EXAMPLES_TYPES_HPP

#include <fwdpp/ts/definitions.hpp>
#include <fwdpp/diploid_population.hpp>
#include <fwdpp/GSLrng_t.hpp>
#include <fwdpp/types/mutation.hpp>
#include <fwdpp/ts/table_collection.hpp>

using ts_examples_poptype = fwdpp::diploid_population<fwdpp::mutation>;
using GSLrng = fwdpp::GSLrng_mt;

struct diploid_metadata
{
    std::size_t individual;
    double time, fitness;
    fwdpp::ts::table_collection::id_type n1, n2;
    diploid_metadata(std::size_t i, double t, double w,
                     fwdpp::ts::table_collection::id_type a,
                     fwdpp::ts::table_collection::id_type b)
        : individual(i), time(t), fitness(w), n1(a), n2(b)
    {
    }
};

#endif
