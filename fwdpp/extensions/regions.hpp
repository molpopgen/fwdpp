///! \file regions.hpp
#ifndef __FWDPP_EXTENSIONS_REGIONS_HPP__
#define __FWDPP_EXTENSIONS_REGIONS_HPP__

#include <cassert>
#include <limits>
#include <stdexcept>
#include <algorithm>
#include <memory>
#include <gsl/gsl_randist.h>
#include <fwdpp/type_traits.hpp>
#include <fwdpp/gsl_discrete.hpp>
#include <fwdpp/simfunctions/recycling.hpp>
#include <fwdpp/extensions/callbacks.hpp>

namespace fwdpp
{
    namespace extensions
    {
        template <typename MutationContainerType, typename GenomeContainerType = void,
                  typename MutationModelSignatureType = typename std::conditional<
                      std::is_void<GenomeContainerType>::value,
                      typename fwdpp::traits::mutation_model<MutationContainerType>,
                      typename fwdpp::traits::mutation_model_haploid_genome<
                          MutationContainerType, GenomeContainerType>>::type>
        class discrete_mut_model
        /*!
         *  A container of mutation models + weights.
         *  When called, specific model functions are applied according
         *  to their weights.
         *
         *  Can be used to model regional variation in mutational
         *  processes and/or mixtures of different distributions
         *  of effect sizes.
         *
         *  fwdpp 0.6.0 generalized the implementation to be
         *  in terms of std::function.
         *
         *  See extensions_regionsTest.cc and
         *  K_linked_regions_extensions.cc for examples.
         * 
         *  \example K_linked_regions_extensions.cc
         */
        {
            static_assert(fwdpp::traits::is_mutation<
                              typename MutationContainerType::value_type>::value,
                          "MutationContainerType::value_type must be a mutation type");

          private:
            std::vector<MutationModelSignatureType> functions;
            std::vector<double> weights;
            fwdpp::gsl_ran_discrete_t_ptr lookup;

          public:
            using function_type = MutationModelSignatureType;
            using mutation_container = MutationContainerType;
            using haploid_genome_container = GenomeContainerType;

            template <typename FunctionContainerType, typename WeightContainerType>
            discrete_mut_model(FunctionContainerType &&functions_,
                               WeightContainerType &&weights_)
                : functions{std::forward<FunctionContainerType>(functions_)},
                  weights{std::forward<WeightContainerType>(weights_)}, lookup{nullptr}
            {
                if (functions.size() != weights.size())
                    {
                        throw std::invalid_argument("number of functions must "
                                                    "equal number of weights");
                    }
                if (!weights.empty())
                    {
                        lookup = fwdpp::gsl_ran_discrete_t_ptr(
                            gsl_ran_discrete_preproc(weights.size(), weights.data()));
                    }
                else
                    {
                        throw std::invalid_argument("weights cannot be empty");
                    }
            }

            discrete_mut_model(const discrete_mut_model &dmm)
                : functions{dmm.functions}, weights{dmm.weights}, lookup{nullptr}
            {
                if (!weights.empty())
                    {
                        lookup = fwdpp::gsl_ran_discrete_t_ptr(
                            gsl_ran_discrete_preproc(weights.size(), weights.data()));
                    }
            }

            template <typename... function_type_args>
            inline typename function_type::result_type
            operator()(const gsl_rng *r, function_type_args &&...args) const
            {
                auto region = gsl_ran_discrete(r, lookup.get());
                return functions.at(region)(std::forward<function_type_args>(args)...);
            }

            const std::vector<double> &
            get_weights() const
            /// Returns const reference to the weights.
            {
                return weights;
            }
        };

        template <typename dmm_type>
        inline typename dmm_type::function_type
        bind_dmm_wrapper(const gsl_rng *r, const dmm_type &dmm, std::true_type)
        {
            return [r, &dmm](flagged_mutation_queue &recycling_bin,
                             typename dmm_type::mutation_container &mutations) {
                return dmm(r, recycling_bin, mutations);
            };
        }

        template <typename dmm_type>
        inline typename dmm_type::function_type
        bind_dmm_wrapper(const gsl_rng *r, const dmm_type &dmm, std::false_type)
        {
            return
                [r, &dmm](const typename dmm_type::haploid_genome_container::value_type
                              &haploid_genome,
                          flagged_mutation_queue &recycling_bin,
                          typename dmm_type::mutation_container &mutations) {
                    return dmm(r, haploid_genome, recycling_bin, mutations);
                };
        }

        /*!
          Convenience function to return a function object
          bound to fwdpp::extensions::discrete_mut_model::operator()

          See extensions_regionsTest.cc and
          K_linked_regions_extensions.cc for examples.
        */
        template <typename dmm_type>
        inline typename dmm_type::function_type
        bind_dmm(const gsl_rng *r, const dmm_type &dmm)
        {
            return bind_dmm_wrapper(
                r, dmm,
                typename std::is_void<
                    typename dmm_type::haploid_genome_container>::type());
        }

        /*! Return a vector of callables bound
         *  to fwdpp::extensions::discrete_mut_model::operator()
         *  See extensions_regionsTest.cc and
         *  K_linked_regions_extensions.cc for examples.
         */
        template <typename dmm_type>
        inline std::vector<typename dmm_type::function_type>
        bind_vec_dmm(const gsl_rng *r, const std::vector<dmm_type> &vdm)
        {
            std::vector<typename dmm_type::function_type> rv;
            for (auto &i : vdm)
                {
                    rv.emplace_back(bind_dmm(r, i));
                }
            return rv;
        }

        class discrete_rec_model
        /*!
          Class allowing the simulation of discrete variation
          in recombination rates along a region.

          \example K_linked_regions_extensions.cc
        */
        {
          private:
            double recrate;
            std::vector<std::function<void(std::vector<double> &)>> functions;
            std::vector<double> weights;
            fwdpp::gsl_ran_discrete_t_ptr lookup;

          public:
            using result_type = std::vector<double>;
            using function_type = std::function<void(std::vector<double> &)>;
            /*!
             */
            template <typename FunctionContainerType, typename WeightContainerType>
            discrete_rec_model(const double recrate_, FunctionContainerType &&functions_,
                               WeightContainerType &&weights_)
                : recrate{recrate_},
                  functions(std::forward<FunctionContainerType>(functions_)),
                  weights(std::forward<WeightContainerType>(weights_)), lookup(nullptr)
            {
                if (functions.size() != weights.size())
                    {
                        throw std::invalid_argument("number of functions must "
                                                    "equal number of weights");
                    }
                if (!weights.empty())
                    {
                        lookup = fwdpp::gsl_ran_discrete_t_ptr(
                            gsl_ran_discrete_preproc(weights.size(), weights.data()));
                    }
                else
                    {
                        throw std::invalid_argument("weights cannot be empty");
                    }
            }

            discrete_rec_model(const discrete_rec_model &drm)
                : recrate{drm.recrate}, functions{drm.functions}, weights{drm.weights},
                  lookup{nullptr}

            {
                if (!weights.empty())
                    {
                        lookup = fwdpp::gsl_ran_discrete_t_ptr(
                            gsl_ran_discrete_preproc(weights.size(), weights.data()));
                    }
            }

            /*!
             */
            inline result_type
            operator()(const gsl_rng *r) const
            {
                unsigned nbreaks = gsl_ran_poisson(r, recrate);
                if (!nbreaks)
                    return {};

                std::vector<double> rv;
                for (unsigned i = 0; i < nbreaks; ++i)
                    {
                        auto region = gsl_ran_discrete(r, lookup.get());
                        assert(region < functions.size());
                        functions.at(region)(rv);
                    }
                std::sort(rv.begin(), rv.end());
                rv.push_back(std::numeric_limits<double>::max());
                return rv;
            }
            const std::vector<double> &
            get_weights() const
            /// Returns const reference to the weights.
            {
                return weights;
            }

            double
            get_recrate() const
            /// Returns the recombination rate.
            {
                return recrate;
            }
        };
    } // namespace extensions
} // namespace fwdpp
#endif
