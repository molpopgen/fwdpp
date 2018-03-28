///! \file regions.hpp
#ifndef __FWDPP_EXTENSIONS_REGIONS_HPP__
#define __FWDPP_EXTENSIONS_REGIONS_HPP__

#include <limits>
#include <cassert>
#include <stdexcept>
#include <algorithm>
#include <memory>
#include <gsl/gsl_randist.h>
#include <fwdpp/type_traits.hpp>
#include <fwdpp/sugar/infsites.hpp>
#include <fwdpp/internal/gsl_discrete.hpp>
#include <fwdpp/internal/recycling.hpp>
#include <fwdpp/extensions/callbacks.hpp>

namespace fwdpp
{
    namespace extensions
    {
        template <
            typename mcont_t, typename gcont_t = void,
            typename mutation_model_signature = typename std::conditional<
                std::is_void<gcont_t>::value,
                typename fwdpp::traits::mutation_model<mcont_t>,
                typename fwdpp::traits::mutation_model_gamete<mcont_t,
                                                              gcont_t>>::type>
        class discrete_mut_model
        /*!
         */
        {
            static_assert(fwdpp::traits::is_mutation<
                              typename mcont_t::value_type>::value,
                          "mcont_t::value_type must be a mutation type");

          private:
            std::vector<mutation_model_signature> functions;
            std::vector<double> weights;
            fwdpp::fwdpp_internal::gsl_ran_discrete_t_ptr lookup;

          public:
            using function_type = mutation_model_signature;
            using mutation_container = mcont_t;
            using gamete_container = gcont_t;
            using recycling_bin_type =
                typename fwdpp::traits::recycling_bin_t<mcont_t>;

            template <typename function_container, typename weight_container>
            discrete_mut_model(function_container &&functions_,
                               weight_container &&weights_)
                : functions{ std::forward<function_container>(functions_) },
                  weights{ std::forward<weight_container>(weights_) }, lookup{
                      nullptr
                  }
            {
                if (functions.size() != weights.size())
                    {
                        throw std::invalid_argument("number of functions must "
                                                    "equal number of weights");
                    }
                if (!weights.empty())
                    {
                        lookup = fwdpp::fwdpp_internal::gsl_ran_discrete_t_ptr(
                            gsl_ran_discrete_preproc(weights.size(),
                                                     weights.data()));
                    }
                else
                    {
                        throw std::invalid_argument("weights cannot be empty");
                    }
            }

            template <typename... function_type_args>
            inline typename function_type::result_type
            operator()(const gsl_rng *r, function_type_args &&... args) const
            {
                auto region = gsl_ran_discrete(r, lookup.get());
                return functions.at(region)(
                    std::forward<function_type_args>(args)...);
            }
        };

        template <typename dmm_type>
        inline typename dmm_type::function_type
        bind_dmm_wrapper(const gsl_rng *r, dmm_type &dmm, std::true_type)
        {
            return
                [r, &dmm](typename dmm_type::recycling_bin_type &recycling_bin,
                          typename dmm_type::mutation_container &mutations) {
                    return dmm(r, recycling_bin, mutations);
                };
        }

        template <typename dmm_type>
        inline typename dmm_type::function_type
        bind_dmm_wrapper(const gsl_rng *r, dmm_type &dmm, std::false_type)
        {
            return
                [r, &dmm](typename dmm_type::recycling_bin_type &recycling_bin,
                          typename dmm_type::gcont_t::value_type &gamete,
                          typename dmm_type::mutation_container &mutations) {
                    return dmm(r, recycling_bin, gamete, mutations);
                };
        }

        /*!
          Convenience function to return a function object
          bound to fwdpp::extensions::discrete_mut_model::operator()

          This simplifies life a lot!

        */
        template <typename dmm_type>
        inline typename dmm_type::function_type
        bind_dmm(const gsl_rng *r, dmm_type &dmm)
        {
            return bind_dmm_wrapper(
                r, dmm,
                typename std::is_void<
                    typename dmm_type::gamete_container>::type());
        }

        /*! Return a vector of callables bound
         *  to fwdpp::extensions::discrete_mut_model::operator()
         */
        template <typename dmm_type>
        inline std::vector<typename dmm_type::function_type>
        bind_vec_dmm(const gsl_rng *r, std::vector<dmm_type> &vdm)
        {
            std::vector<typename dmm_type::function_type> rv;
            for (auto &i : vdm)
                {
                    rv.emplace_back(bind_dmm(r, i));
                }
            return rv;
        }

        struct discrete_rec_model_data
        {
            const gsl_rng *r;
            const double recrate;
            std::vector<double> beg, end, weight;
            discrete_rec_model_data(const gsl_rng *r_, const double recrate_,
                                    const std::vector<double> &&b,
                                    const std::vector<double> &&e,
                                    const std::vector<double> &&w)
                : r{ r_ }, recrate{ recrate_ }, beg(std::move(b)),
                  end(std::move(e)), weight(std::move(w))
            {
                if (beg.size() != end.size() || beg.size() != weight.size())
                    {
                        throw std::invalid_argument(
                            "input vectors must all be the same size");
                    }
            }
        };

        class discrete_rec_model
        /*!
          Class allowing the simulation of discrete variation
          in recombination rates along a region.
        */
        {
          private:
            using result_type = std::vector<double>;
            std::unique_ptr<discrete_rec_model_data> data;
            fwdpp::fwdpp_internal::gsl_ran_discrete_t_ptr lookup;
            void
            assign_weights()
            {
                if (data->weight.size())
                    lookup = fwdpp::fwdpp_internal::gsl_ran_discrete_t_ptr(
                        gsl_ran_discrete_preproc(data->weight.size(),
                                                 data->weight.data()));
            }

          public:
            /*!
                          \param r A random number generator
                          \param recrate Expected number of breakpoints per
              diploid
              \param __beg Region beginnings
              \param __end Region ends
              \param __weight Region weights
            */
            discrete_rec_model(const gsl_rng *r, const double recrate,
                               const std::vector<double> __beg,
                               const std::vector<double> __end,
                               const std::vector<double> __weight)
                : data(new discrete_rec_model_data(
                      r, recrate, std::move(__beg), std::move(__end),
                      std::move(__weight)))
            {
                assign_weights();
            }

            discrete_rec_model(const discrete_rec_model &drm)
                : data(new discrete_rec_model_data(*drm.data))
            {
                assign_weights();
            }

            /*!
              Returns a position from a region that is chosen based on
              region
              weights.

              \precondition If recrate == 0, then beg and end cannot be
              empty.
              It is up to the calling environment to make sure
              that this cannot happen.  This is checked in debug mode via
              the
              assert macro.
            */
            inline result_type
            operator()() const
            {
                assert(!(data->recrate == 0.
                         && (data->beg.empty() || data->end.empty())));
                auto nbreaks = gsl_ran_poisson(data->r, data->recrate);
                if (!nbreaks)
                    return {};

                result_type rv;
                for (unsigned i = 0; i < nbreaks; ++i)
                    {
                        size_t region
                            = gsl_ran_discrete(data->r, lookup.get());
                        rv.push_back(gsl_ran_flat(data->r, data->beg[region],
                                                  data->end[region]));
                    }
                std::sort(rv.begin(), rv.end());
                rv.push_back(std::numeric_limits<double>::max());
                return rv;
            }
        };
    }
}
#endif
