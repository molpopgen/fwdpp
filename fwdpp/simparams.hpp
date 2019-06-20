#ifndef FWDPP_SIM_PARAMS_HPP
#define FWDPP_SIM_PARAMS_HPP

#include <vector>
#include <queue>
#include <cstdint>
#include <utility>
#include <type_traits>
#include <gsl/gsl_rng.h>
#include <fwdpp/type_traits.hpp>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/simfunctions/recycling.hpp>

namespace fwdpp
{
    struct mendel
    /// Returns 1 one-half of the time, indicating a swapping
    /// of parental haploid_genomes
    /// \version 0.7.4 Added to fwdpp
    {
        mendel(){};
        inline int
        operator()(const gsl_rng* r, std::size_t /*parental_haploid_genome_1*/,
                   std::size_t /*parental_haploid_genome_2*/) const
        {
            return gsl_rng_uniform(r) < 0.5;
        }
    };

    template <typename genetic_value, typename mutation_function,
              typename recombination_function,
              typename interlocus_recombination_function,
              typename parent_haploid_genome_swapping_function>
    struct genetic_parameters
    /// \brief Hold types for genetic operations
    ///
    /// This is a simplistic struct to encapsulate the
    /// types needed for "doing the genetics" in a simulation.
    ///
    /// There is no effort to ensure that the types reflect
    /// valid policies upon instantiation or construction.
    ///
    /// For a single-locus simulation, interlocus_recombination_function
    /// is expected to be std::nullptr_t
    ///
    /// The heavy use of lambdas in fwdpp means that it is, in practice,
    /// nearly impossible to know the template parameter types.  Thus,
    /// objects of this type are made by calling fwdpp::make_genetic_parameters.
    ///
    /// \version 0.7.4 Added to library
    {
        const genetic_value gvalue;
        const mutation_function generate_mutations;
        const recombination_function generate_breakpoints;
        const interlocus_recombination_function interlocus_recombination;
        const parent_haploid_genome_swapping_function haploid_genome_swapper;
        flagged_mutation_queue mutation_recycling_bin;
        flagged_haploid_genome_queue haploid_genome_recycling_bin;
        std::vector<uint_t> neutral;
        std::vector<uint_t> selected;

        template <typename gv, typename mut, typename rec, typename swapper>
        genetic_parameters(gv&& gvalue_param, mut&& generate_mutations_param,
                           rec&& generate_breakpoints_param,
                           swapper&& haploid_genome_swapper_param)
            : gvalue{ std::forward<gv>(gvalue_param) },
              generate_mutations{ std::forward<mut>(
                  generate_mutations_param) },
              generate_breakpoints{ std::forward<rec>(
                  generate_breakpoints_param) },
              interlocus_recombination{},
              haploid_genome_swapper{ std::forward<swapper>(haploid_genome_swapper_param) },
              mutation_recycling_bin{ empty_mutation_queue() },
              haploid_genome_recycling_bin{ empty_haploid_genome_queue() }, neutral{},
              selected{}
        {
        }

        template <typename gv, typename mut, typename rec, typename irec,
                  typename swapper>
        genetic_parameters(gv&& gvalue_param, mut&& generate_mutations_param,
                           rec&& generate_breakpoints_param,
                           irec&& interlocus_recombination_param,
                           swapper&& haploid_genome_swapper_param)
            : gvalue{ std::forward<gv>(gvalue_param) },
              generate_mutations{ std::forward<mut>(
                  generate_mutations_param) },
              generate_breakpoints{ std::forward<rec>(
                  generate_breakpoints_param) },
              interlocus_recombination{ std::forward<irec>(
                  interlocus_recombination_param) },
              haploid_genome_swapper{ std::forward<swapper>(haploid_genome_swapper_param) },
              mutation_recycling_bin{ empty_mutation_queue() },
              haploid_genome_recycling_bin{ empty_haploid_genome_queue() }, neutral{},
              selected{}
        {
        }
    };

    template <typename genetic_value, typename mutation_function,
              typename recombination_function>
    inline genetic_parameters<
        typename std::remove_reference<genetic_value>::type,
        typename std::remove_reference<mutation_function>::type,
        typename std::remove_reference<recombination_function>::type,
        std::nullptr_t, mendel>
    make_genetic_parameters(
        genetic_value&& gvalue_param,
        mutation_function&& generate_mutations_param,
        recombination_function&& generate_breakpoints_param)
    /// Create a fwdpp::genetic_parameters for a single-locus simulation.
    {
        return genetic_parameters<
            typename std::remove_reference<genetic_value>::type,
            typename std::remove_reference<mutation_function>::type,
            typename std::remove_reference<recombination_function>::type,
            std::nullptr_t, mendel>(
            std::forward<genetic_value>(gvalue_param),
            std::forward<mutation_function>(generate_mutations_param),
            std::forward<recombination_function>(generate_breakpoints_param),
            mendel());
    }

    template <typename genetic_value, typename mutation_function,
              typename recombination_function, typename swapper>
    inline genetic_parameters<
        typename std::remove_reference<genetic_value>::type,
        typename std::remove_reference<mutation_function>::type,
        typename std::remove_reference<recombination_function>::type,
        std::nullptr_t, typename std::remove_reference<swapper>::type>
    make_genetic_parameters_with_swapper(
        genetic_value&& gvalue_param,
        mutation_function&& generate_mutations_param,
        recombination_function&& generate_breakpoints_param,
        swapper&& swapper_param)
    /// Create a fwdpp::genetic_parameters for a single-locus simulation.
    {
        return genetic_parameters<
            typename std::remove_reference<genetic_value>::type,
            typename std::remove_reference<mutation_function>::type,
            typename std::remove_reference<recombination_function>::type,
            std::nullptr_t, typename std::remove_reference<swapper>::type>(
            std::forward<genetic_value>(gvalue_param),
            std::forward<mutation_function>(generate_mutations_param),
            std::forward<recombination_function>(generate_breakpoints_param),
            std::forward<swapper>(swapper_param));
    }
} // namespace fwdpp

#endif

