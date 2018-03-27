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
        template <typename mcont_t,
                  typename mutation_model_signature =
                      typename fwdpp::traits::mutation_model<mcont_t>>
        class discrete_mut_model
        /*!
        */
        {
          private:
            std::vector<mutation_model_signature> functions;
            std::vector<double> weights;
            fwdpp::fwdpp_internal::gsl_ran_discrete_t_ptr lookup;

          public:
            using function_type = mutation_model_signature;

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

            inline typename function_type::result_type
            operator()(const gsl_rng *r,
                       typename fwdpp::traits::recycling_bin_t<mcont_t>
                           &mutation_recycling_bin,
                       mcont_t &mutations) const
            {
                auto region = gsl_ran_discrete(r, lookup.get());
                return functions[region](mutation_recycling_bin, mutations);
            }
        };

        /*!
          Convenience function to return a function object
          bound to fwdpp::extensions::discrete_mut_model::operator()

          This simplifies life a lot!

          See unit test extensions.cc for an example of use.
        */
        // template <typename mcont_t, typename lookup_t, class... Args>
        // inline traits::mutation_model<mcont_t>
        // bind_dmm(const discrete_mut_model &dm, mcont_t &, lookup_t
        // &mut_lookup,
        //          Args &&... args)
        // {
        //     return std::bind(&discrete_mut_model::
        //                      operator()<traits::recycling_bin_t<mcont_t>,
        //                                 lookup_t, mcont_t>,
        //                      &dm, std::placeholders::_1,
        //                      std::placeholders::_2,
        //                      std::forward<Args>(args)...,
        //                      std::ref(mut_lookup));
        // }

        // /*! Return a vector of callables bount
        //  *  to fwdpp::extensions::discrete_mut_model::operator()
        //  */
        // template <typename mcont_t, typename lookup_t, class... Args>
        // inline std::vector<traits::mutation_model<mcont_t>>
        // bind_vec_dmm(const std::vector<discrete_mut_model> &vdm,
        //              mcont_t &mutations, lookup_t &mut_lookup,
        //              const gsl_rng *r,
        //              const std::vector<double> &neutral_mutrates,
        //              const std::vector<double> &selected_mutrates,
        //              Args &&... args)
        // {
        //     if (vdm.size() != neutral_mutrates.size()
        //         || vdm.size() != selected_mutrates.size())
        //         {
        //             throw std::invalid_argument(
        //                 "container sizes must all be equal");
        //         }
        //     std::vector<decltype(
        //         bind_dmm(vdm[0], mutations, mut_lookup, r,
        //         neutral_mutrates[0],
        //                  selected_mutrates[0],
        //                  std::forward<Args>(args)...))>
        //         rv;
        //     std::size_t i = 0;
        //     for (auto &&dm : vdm)
        //         {
        //             rv.emplace_back(bind_dmm(
        //                 dm, mutations, mut_lookup, r, neutral_mutrates[i],
        //                 selected_mutrates[i], std::forward<Args>(args)...));
        //             ++i;
        //         }
        //     return rv;
        // }

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
