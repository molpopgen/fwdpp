#ifndef TS_EXAMPLES_COMMON_HPP__
#define TS_EXAMPLES_COMMON_HPP__

#include <cstdint>
#include <string>
#include <functional>
#include <boost/program_options.hpp>
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

struct options
{
    fwdpp::uint_t N, gcint;
    double theta, rho, mean, shape, mu, scoeff, dominance, scaling;
    unsigned seed;
    int ancient_sampling_interval, ancient_sample_size, nsam, nloci;
    bool leaf_test, matrix_test, preserve_fixations;
    std::string filename, sfsfilename;
    options();
};

boost::program_options::options_description generate_main_options(options &o);
boost::program_options::options_description generate_dfe_options(options &o);
boost::program_options::options_description
generate_testing_options(options &);
void validate_primary_options(const options &);

std::function<double()> make_dfe(const fwdpp::uint_t N,
                                 const fwdpp::GSLrng_mt &r, const double mean,
                                 const double shape, const double scoeff);

void matrix_runtime_test(const fwdpp::ts::table_collection &tables,
                         const std::vector<fwdpp::ts::TS_NODE_INT> &samples,
                         const std::vector<fwdpp::popgenmut> &mutations,
                         const std::vector<fwdpp::uint_t> &mcounts);

void
expensive_leaf_test(const fwdpp::ts::table_collection &tables,
                    const std::vector<fwdpp::ts::TS_NODE_INT> &sample_list);

void test_serialization(const fwdpp::ts::table_collection &tables,
                        const std::string &filename);
#endif
