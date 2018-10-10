/*! \include wfts.cet
 */

#include <cstdio>
#include <iostream>
#include <cassert>
#include <fwdpp/ts/table_collection.hpp>
#include <fwdpp/ts/table_simplifier.hpp>
#include <fwdpp/sugar/popgenmut.hpp>
#include <fwdpp/sugar/slocuspop.hpp>
#include <boost/program_options.hpp>

#include "evolve_generation_ts.hpp"

namespace po = boost::program_options;
using poptype = fwdpp::slocuspop<fwdpp::popgenmut>;

int
main(int argc, char** argv)
{
    fwdpp::uint_t N, gcint;
    double theta,rho,mean,shape,mu;
    po::options_description options("Usage");
    options.add_options()("help", "Display help")
        ("N", po::value<unsigned>(&N), "Diploid population size")
        ("gc", po::value<unsigned>(&gcint),"Simplification interval")
        ("theta", po::value<double>(&theta), "4Nu")
        ("rho", po::value<double>(&rho), "4Nr")
        ("mu", po::value<double>(&mu), "mutation rate to selected variants")
        ("mean", po::value<double>(&mean),
        "Mean of Gamma distribution of selection coefficients")
        ("shape", po::value<double>(&shape),
        "Shape of Gamma distribution of selection coefficients");
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, options), vm);
    po::notify(vm);

    if (vm.count("help"))
        {
            std::cout << options << '\n';
            std::exit(1);
        }
    fwdpp::ts::table_collection tables(1.0);
    fwdpp::ts::table_simplifier(1.0);
}
