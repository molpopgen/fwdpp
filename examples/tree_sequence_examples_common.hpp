#ifndef TS_EXAMPLES_COMMON_HPP__
#define TS_EXAMPLES_COMMON_HPP__

#include <cstdint>
#include <string>
#include <functional>
#include <boost/program_options.hpp>
#include <fwdpp/simfunctions/recycling.hpp>
#include <fwdpp/ts/std_table_collection.hpp>
#include "tree_sequence_examples_types.hpp"

struct options
{
    fwdpp::uint_t N, gcint;
    double theta, rho, mean, shape, mu, scoeff, dominance, scaling;
    unsigned seed;
    int ancient_sampling_interval, ancient_sample_size, nsam;
    bool leaf_test, matrix_test, visit_sites_test, preserve_fixations;
    std::string filename, sfsfilename;
    options();
};

int apply_neutral_mutations(const options &o, const GSLrng &rng,
                            fwdpp::ts::std_table_collection &tables,
                            ts_examples_poptype &pop,
                            fwdpp::flagged_mutation_queue &mutation_recycling_bin);

boost::program_options::options_description generate_main_options(options &o);
boost::program_options::options_description generate_dfe_options(options &o);
boost::program_options::options_description generate_testing_options(options &);
void validate_primary_options(const options &);

std::function<double()> make_dfe(const fwdpp::uint_t N, const GSLrng &r,
                                 const double mean, const double shape,
                                 const double scoeff);

void execute_matrix_test(const options &, const ts_examples_poptype &,
                         const fwdpp::ts::std_table_collection &,
                         const std::vector<fwdpp::ts::table_index_t> &,
                         const std::vector<fwdpp::ts::table_index_t> &);

void execute_expensive_leaf_test(
    const options &o, const fwdpp::ts::std_table_collection &tables,
    const std::vector<fwdpp::ts::table_index_t> &samples,
    const std::vector<fwdpp::ts::table_index_t> &preserved_nodes);

void execute_serialization_test(const options &,
                                const fwdpp::ts::std_table_collection &);

void test_serialization(const fwdpp::ts::std_table_collection &tables,
                        const std::string &filename);

void visit_sites_test(const options &o, const ts_examples_poptype &pop,
                      const fwdpp::ts::std_table_collection &tables,
                      const std::vector<fwdpp::ts::table_index_t> &samples,
                      const std::vector<fwdpp::ts::table_index_t> &preserved_nodes);

void write_sfs(const options &o, const GSLrng &rng,
               const fwdpp::ts::std_table_collection &tables,
               const std::vector<fwdpp::ts::table_index_t> &samples);

#endif
