#include <vector>
#include <memory>
#include <fwdpp/GSLrng_t.hpp>
#include <fwdpp/recbinder.hpp>
#include <fwdpp/genetic_map/poisson_interval.hpp>
#include <fwdpp/genetic_map/binomial_point.hpp>
#include <fwdpp/genetic_map/genetic_map.hpp>

/* Below, we have two methods to create a fwdpp::genetic_map.
 * Both create the same map.  Each map has 3 "regions":
 * 1. A region that generates a Poisson number of breakpoints on the continuous interval [0, 10)
 * 2. A region that generates a breakpoint at position 10 with probability 0.5.
 * 3. A region that generates a Poisson number of breakpoints on the continuous interval [10, 20)
 *
 * This model is therefore two regions separated by 50 centiMorgans.
 *
 * methods 1 and 2 are the more efficient of the three, requiring fewer temporary objects.
 * method 3 is useful for situations where objects must be accessed directly.  For example,
 * when constructing Python interfaces via pybind11, one cannot access a unique_ptr.
 * It is unlikely that any method would be a measurable performance bottleneck.
 */

fwdpp::genetic_map
method1()
{
    std::vector<std::unique_ptr<fwdpp::genetic_map_unit>> callbacks;
    callbacks.emplace_back(std::make_unique<fwdpp::poisson_interval>(0, 10, 1e-3));
    callbacks.emplace_back(std::make_unique<fwdpp::binomial_point>(10., 0.5));
    callbacks.emplace_back(std::make_unique<fwdpp::poisson_interval>(10., 20, 1e-3));
    return fwdpp::genetic_map(std::move(callbacks));
}

fwdpp::genetic_map
method2()
{
    fwdpp::genetic_map gmap;
    gmap.add_callback(std::make_unique<fwdpp::poisson_interval>(0, 10, 1e-3));
    gmap.add_callback(std::make_unique<fwdpp::binomial_point>(10., 0.5));
    gmap.add_callback(std::make_unique<fwdpp::poisson_interval>(10., 20, 1e-3));
    return gmap;
}

fwdpp::genetic_map
method3()
{
    fwdpp::genetic_map gmap;
    gmap.add_callback(fwdpp::poisson_interval(0, 10, 1e-3));
    gmap.add_callback(fwdpp::binomial_point(10., 0.5));
    gmap.add_callback(fwdpp::poisson_interval(10., 20, 1e-3));
    return gmap;
}

int
main(int /*argc*/, char** /*argv*/)
{
    auto map_method1 = method1();
    auto map_method2 = method2();
    auto map_method3 = method3();

    fwdpp::GSLrng_mt rng(42);

    // A simulation expects a genetic map to be equivalent to
    // std::function<std::vector<double>(void)>.
    // Our variables map_method1 and map_method2 are not.
    // The function fwdpp::recbinder rebinds them so that
    // they conform to API requirements.

    auto bound_map1 = fwdpp::recbinder(std::cref(map_method1), rng.get());
    auto bound_map2 = fwdpp::recbinder(std::cref(map_method2), rng.get());
    auto bound_map3 = fwdpp::recbinder(std::cref(map_method3), rng.get());

    // Generate breakpoints from our model.
    // You'd use "auto" here, but we write the exact return
    // type for documentation purposes.
    std::vector<double> breakpoints = bound_map1();
    breakpoints = bound_map2();
    breakpoints = bound_map3();
}

