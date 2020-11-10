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

    template <typename GeneticValueType, typename MutationFunctionType,
              typename RecombinationFunctionType, typename SwappingFunctionType>
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
    /// \version 0.9.0 gvalue no longer const
    {
        using genetic_value = GeneticValueType;
        using mutation_function = MutationFunctionType;
        using recombination_function = RecombinationFunctionType;
        using swapping_function = SwappingFunctionType;
        genetic_value gvalue;
        const mutation_function generate_mutations;
        const recombination_function generate_breakpoints;
        const swapping_function haploid_genome_swapper;
        flagged_mutation_queue mutation_recycling_bin;
        flagged_haploid_genome_queue haploid_genome_recycling_bin;
        std::vector<uint_t> neutral;
        std::vector<uint_t> selected;

        template <typename gv, typename mut, typename rec, typename swapper>
        genetic_parameters(gv&& gvalue_param, mut&& generate_mutations_param,
                           rec&& generate_breakpoints_param,
                           swapper&& haploid_genome_swapper_param)
            : gvalue{std::forward<gv>(gvalue_param)},
              generate_mutations{std::forward<mut>(generate_mutations_param)},
              generate_breakpoints{std::forward<rec>(generate_breakpoints_param)},
              haploid_genome_swapper{
                  std::forward<swapper>(haploid_genome_swapper_param)},
              mutation_recycling_bin{empty_mutation_queue()},
              haploid_genome_recycling_bin{empty_haploid_genome_queue()}, neutral{},
              selected{}
        {
        }
    };

    template <typename GeneticValueType, typename MutationFunctionType,
              typename RecombinationFunctionType, typename SwappingFunctionType>
    inline genetic_parameters<
        typename std::remove_reference<GeneticValueType>::type,
        typename std::remove_reference<MutationFunctionType>::type,
        typename std::remove_reference<RecombinationFunctionType>::type,
        typename std::remove_reference<SwappingFunctionType>::type>
    make_genetic_parameters(GeneticValueType&& gvalue_param,
                            MutationFunctionType&& generate_mutations_param,
                            RecombinationFunctionType&& generate_breakpoints_param,
                            SwappingFunctionType&& swapper_param)
    /// Create a fwdpp::genetic_parameters using a "custom swapper
    /// instead of fwdpp::mendel.
    {
        return genetic_parameters<
            typename std::remove_reference<GeneticValueType>::type,
            typename std::remove_reference<MutationFunctionType>::type,
            typename std::remove_reference<RecombinationFunctionType>::type,
            typename std::remove_reference<SwappingFunctionType>::type>(
            std::forward<GeneticValueType>(gvalue_param),
            std::forward<MutationFunctionType>(generate_mutations_param),
            std::forward<RecombinationFunctionType>(generate_breakpoints_param),
            std::forward<SwappingFunctionType>(swapper_param));
    }

    template <typename GeneticValueType, typename MutationFunctionType,
              typename RecombinationFunctionType>
    inline genetic_parameters<
        typename std::remove_reference<GeneticValueType>::type,
        typename std::remove_reference<MutationFunctionType>::type,
        typename std::remove_reference<RecombinationFunctionType>::type, mendel>
    make_genetic_parameters(GeneticValueType&& gvalue_param,
                            MutationFunctionType&& generate_mutations_param,
                            RecombinationFunctionType&& generate_breakpoints_param)
    /// Create a fwdpp::genetic_parameters using fwdpp::mendel
    {
        return make_genetic_parameters(
            std::forward<GeneticValueType>(gvalue_param),
            std::forward<MutationFunctionType>(generate_mutations_param),
            std::forward<RecombinationFunctionType>(generate_breakpoints_param),
            mendel());
    }
} // namespace fwdpp

#endif

