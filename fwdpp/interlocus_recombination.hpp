//
// Copyright (C) 2017 Kevin Thornton <krthornt@uci.edu>
//
// This file is part of fwdpp.
//
// fwdpp is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// fwdpp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with fwdpp.  If not, see <http://www.gnu.org/licenses/>.
//

/*! \file interlocus_recombination.hpp
 Defines convenience functions to aid multi-locus/region simulations
*/

#ifndef FWDPP_INTERLOCUS_RECOMBINATION_HPP__
#define FWDPP_INTERLOCUS_RECOMBINATION_HPP__

#include <cstdint>
#include <vector>
#include <functional>
#include <fwdpp/fwd_functional.hpp>

namespace KTfwd
{
    inline std::vector<std::function<unsigned(void)>>
    make_poisson_interlocus_rec(const gsl_rng* r, const double* means,
                                const std::size_t n)
    /// \brief Create a vector of callbacks bound to model interlocus recombination
    ///
    /// Example use for two loci with recombination occuring
    /// at rate 1e-3 between them:
    ///
    /// \code
    /// std::vector<double> recrates_bw_loci{1e-3};
    /// auto interlocus_rec =
    /// KTfwd::make_poisson_interlocus_rec(r,recrates_bw_loci.data(),recrates_bw_loci.size());
    /// \endcode
    ///
    /// \ingroup mlocus
    {
        std::vector<std::function<unsigned(void)>> rv;
        for (std::size_t i = 0; i < n; ++i)
            {
                rv.emplace_back(std::bind(gsl_ran_poisson, r, means[i]));
            }
        return rv;
    }

    inline std::vector<std::function<unsigned(void)>>
    make_binomial_interlocus_rec(const gsl_rng* r, const double* distances,
                                 const std::size_t n)
    /// \brief Create a vector of callbacks bound to model interlocus recombination
    ///
    /// Example of a three locus system where loci 0 and 1
    /// are 25cM apart and locus 2 is unlinked (50cM) from
    /// 0 and 1:
    ///
    /// \code
    /// std::vector<double> recrates{0.25,0.5};
    /// auto interlocus_rec =
    /// KTfwd::make_binomial_interlocus_rec(r,recrates.data(),recrates.size());
    /// \endcode
    ///
    /// \ingroup mlocus
    {
        std::vector<std::function<unsigned(void)>> rv;
        for (std::size_t i = 0; i < n; ++i)
            {
                rv.emplace_back(
                    std::bind(gsl_ran_binomial, r, distances[i], 1));
            }
        return rv;
    }
}

#endif
