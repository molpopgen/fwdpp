/*! \include wfts.cet
 */

#include <cstdio>
#include <iostream>
#include <cassert>
#include <fwdpp/ts/table_collection.hpp>
#include <fwdpp/ts/table_simplifier.hpp>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/sugar/popgenmut.hpp>
#include <fwdpp/sugar/slocuspop.hpp>
#include <fwdpp/extensions/callbacks.hpp>
#include <fwdpp/poisson_xover.hpp>
#include <fwdpp/recbinder.hpp>
#include <boost/program_options.hpp>

#include "evolve_generation_ts.hpp"

namespace po = boost::program_options;
using poptype = fwdpp::slocuspop<fwdpp::popgenmut>;
using GSLrng = fwdpp::GSLrng_t<fwdpp::GSL_RNG_MT19937>;

int
main(int argc, char **argv)
{
    fwdpp::uint_t N, gcint = 100;
    double theta, rho, mean, shape, mu;
    unsigned seed = 42;
    po::options_description options("Usage");
    options.add_options()("help", "Display help")("N", po::value<unsigned>(&N),
                                                  "Diploid population size")(
        "gc", po::value<unsigned>(&gcint),
        "Simplification interval. Default is 100 generations.")(
        "theta", po::value<double>(&theta),
        "4Nu")("rho", po::value<double>(&rho), "4Nr")(
        "mu", po::value<double>(&mu), "mutation rate to selected variants")(
        "mean", po::value<double>(&mean),
        "Mean of Gamma distribution of selection coefficients")(
        "shape", po::value<double>(&shape),
        "Shape of Gamma distribution of selection coefficients")(
        "seed", po::value<unsigned>(&seed),
        "Random number seed. Default is 42");
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, options), vm);
    po::notify(vm);

    if (vm.count("help"))
        {
            std::cout << options << '\n';
            std::exit(1);
        }

    GSLrng rng(seed);

    poptype pop(2 * N);
    fwdpp::ts::table_collection tables(2 * pop.diploids.size(), 0, 0, 1.0);
    fwdpp::ts::table_simplifier(1.0);
    unsigned generation = 1;
    double recrate = rho / static_cast<double>(4 * N);
    const auto recmap
        = fwdpp::recbinder(fwdpp::poisson_xover(recrate, 0., 1.), rng.get());

    fwdpp::extensions::gamma dfe(mean, shape);
    const auto mmodel = [&pop, &rng, &generation,
                         dfe](std::queue<std::size_t> &recbin,
                              poptype::mcont_t &mutations) {
        return fwdpp::infsites_popgenmut(
            recbin, mutations, rng.get(), pop.mut_lookup, generation, 0.0,
            [&rng, dfe]() { return gsl_rng_uniform(rng.get()); },
            [&rng, dfe]() { return dfe(rng.get()); }, []() { return 1.0; });
    };
}
