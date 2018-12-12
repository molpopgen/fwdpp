#ifndef TS_EXAMPLES_COMMON_HPP__
#define TS_EXAMPLES_COMMON_HPP__

#include <cstdint>
#include <string>
#include <functional>
#include <fwdpp/popgenmut.hpp>
#include <fwdpp/GSLrng_t.hpp>
#include <fwdpp/ts/table_collection.hpp>

struct diploid_metadata
{
    std::size_t individual;
    double time, fitness;
    fwdpp::ts::TS_NODE_INT n1, n2;
    diploid_metadata(std::size_t i, double t, double w,
                     fwdpp::ts::TS_NODE_INT a, fwdpp::ts::TS_NODE_INT b)
        : individual(i), time(t), fitness(w), n1(a), n2(b)
    {
    }
};

std::function<double()>
make_dfe(const fwdpp::uint_t N, const fwdpp::GSLrng_mt &r, const double mean,
         const double shape, const double scoeff);


void
matrix_runtime_test(const fwdpp::ts::table_collection &tables,
                    const std::vector<fwdpp::ts::TS_NODE_INT> &samples,
                    const std::vector<fwdpp::popgenmut> &mutations,
                    const std::vector<fwdpp::uint_t> &mcounts);

void
expensive_leaf_test(const fwdpp::ts::table_collection &tables,
                    const std::vector<fwdpp::ts::TS_NODE_INT> &sample_list);

void
test_serialization(const fwdpp::ts::table_collection &tables,
                   const std::string &filename);
#endif
