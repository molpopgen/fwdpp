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

namespace KTfwd
{
    namespace extensions
    {
        struct discrete_mut_model_data
        {
            std::vector<double> nbeg, nend, sbeg, send, nweights, sweights;
            std::vector<shmodel> shmodels;
            using vec_xtra = std::vector<decltype(KTfwd::mutation_base::xtra)>;
            vec_xtra nlabels, slabels;
            discrete_mut_model_data(
                std::vector<double> __nbeg, std::vector<double> __nend,
                std::vector<double> __nweights, // the weights
                std::vector<double> __sbeg, std::vector<double> __send,
                std::vector<double> __sweights, // the weights
                std::vector<shmodel> __shmodels, vec_xtra __nlabels,
                vec_xtra __slabels)
                : nbeg(std::move(__nbeg)), nend(std::move(__nend)),
                  sbeg(std::move(__sbeg)), send(std::move(__send)),
                  nweights(std::move(__nweights)),
                  sweights(std::move(__sweights)),
                  shmodels(std::move(__shmodels)),
                  nlabels(std::move(__nlabels)), slabels(std::move(__slabels))
            {
                if (nbeg.size() != nend.size()
                    || nbeg.size() != nweights.size())
                    {
                        throw std::invalid_argument(
                            "input vectors must all be the same size");
                    }
                if (sbeg.size() != send.size()
                    || sbeg.size() != sweights.size())
                    {
                        throw std::invalid_argument(
                            "input vectors must all be the same size");
                    }
                if (!sweights.empty() && sweights.size() != shmodels.size())
                    {
                        throw std::invalid_argument(
                            "incorrect number of shmodels");
                    }
                if (slabels.empty())
                    {
                        slabels.resize(sbeg.size(), 0);
                    }
                if (nlabels.empty())
                    {
                        nlabels.resize(nend.size(), 0);
                    }
            }
        };

        class discrete_mut_model
        /*!
          Class allowing the simulation of discrete variation
          in mutation models along a region.

          Implemented in terms of KTfwd::popgenmut.

          Intended to be used with the callbacks in
          fwdpp/extensions/callbacks.hpp
        */
        {
          private:
            std::unique_ptr<discrete_mut_model_data> data;
            KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr nlookup, slookup;
            inline void
            assign_weights()
            {
                if (data->nweights.size())
                    nlookup = KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(
                        gsl_ran_discrete_preproc(data->nweights.size(),
                                                 data->nweights.data()));
                if (data->sweights.size())
                    slookup = KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(
                        gsl_ran_discrete_preproc(data->sweights.size(),
                                                 data->sweights.data()));
            }

            /*!
              Helper function.  Not to be called externally.
            */
            template <typename lookup_table_t>
            inline double
            posmaker(const gsl_rng *r, const double &beg, const double &end,
                     lookup_table_t &lookup) const
            {
                double pos = gsl_ran_flat(r, beg, end);
                while (lookup.find(pos) != lookup.end())
                    {
                        pos = gsl_ran_flat(r, beg, end);
                    }
                lookup.insert(pos);
                return pos;
            }

          public:
            using result_type = std::size_t;

            /*!
              \param __nbeg Positions of beginnings of 'neutral' regions
              \param __nend Positions of ends of 'neutral' regions
              \param nweights Weights on 'neutral' regions
              \param __sbeg Positions of beginnings of 'selected' regions
              \param __send Positions of ends of 'selected' regions
              \param sweights Weights on 'selected' regions
              \param __shmodels Vector of KTfwd::experimenta::shmodel
            */
            discrete_mut_model(std::vector<double> __nbeg,
                               std::vector<double> __nend,
                               std::vector<double> nweights, // the weights
                               std::vector<double> __sbeg,
                               std::vector<double> __send,
                               std::vector<double> sweights, // the weights
                               std::vector<shmodel> __shmodels,
                               discrete_mut_model_data::vec_xtra __nlabels
                               = discrete_mut_model_data::vec_xtra(),
                               discrete_mut_model_data::vec_xtra __slabels
                               = discrete_mut_model_data::vec_xtra())
                : data(new discrete_mut_model_data(
                      std::move(__nbeg), std::move(__nend),
                      std::move(nweights), std::move(__sbeg),
                      std::move(__send), std::move(sweights),
                      std::move(__shmodels), std::move(__nlabels),
                      std::move(__slabels)))
            {
                assign_weights();
            }

            discrete_mut_model(const discrete_mut_model &dmm)
                : data(new discrete_mut_model_data(*dmm.data))
            {
                assign_weights();
            }

            /*!
              Return a KTfwd::popgenmut

              \param r A gsl_rng
              \param nmu Neutral mutation rate (per gamete, per generation)
              \param smu Selected mutation rate (per gamete, per
              generation)
              \param generation The current generation
              \param recycling_bin A recycling bin for mutations
              \param mutations A container of mutations
              \param lookup Lookup table for mutations

              \precondition If nmu == 0, then nbeg/nend cannot be empty.
              Similarly,
              if smu == 0, then sbeg,end cannot be empty.  These conditions
              are
              checked in debug
              mode via the assert macro.  It is up to the calling
              environment
              to prevent this situation
              from arising.  Also, nmu+smu must be > 0, and is also checked
              by
              assert.
            */
            template <typename queue_t, typename lookup_table_t,
                      typename mcont_t>
            inline result_type
            make_mut(queue_t &recycling_bin, mcont_t &mutations,
                     const gsl_rng *r, const double nmu, const double smu,
                     const unsigned *generation, lookup_table_t &lookup) const
            {
                assert(nmu + smu > 0.);
                bool is_neutral
                    = (gsl_rng_uniform(r) < nmu / (nmu + smu)) ? true : false;
                if (is_neutral)
                    {
                        assert(!data->nbeg.empty());
                        assert(!data->nend.empty());
                        size_t region = gsl_ran_discrete(r, nlookup.get());
                        double pos = posmaker(r, data->nbeg[region],
                                              data->nend[region], lookup);
                        return fwdpp_internal::recycle_mutation_helper(
                            recycling_bin, mutations, pos, 0., 0., *generation,
                            data->nlabels[region]);
                    }
                assert(!data->sbeg.empty());
                assert(!data->send.empty());
                size_t region = gsl_ran_discrete(r, slookup.get());
                double pos = posmaker(r, data->sbeg[region],
                                      data->send[region], lookup);
                return fwdpp_internal::recycle_mutation_helper(
                    recycling_bin, mutations, pos, data->shmodels[region].s(r),
                    data->shmodels[region].h(r), *generation,
                    data->slabels[region]);
            }
        };

        /*!
          Convenience function to return a function object
          bound to KTfwd::extensions::discrete_mut_model::make_mut

          This simplifies life a lot!

          See unit test extensions.cc for an example of use.
        */
        template <typename mcont_t, typename lookup_t, class... Args>
        inline auto
        bind_dmm(const discrete_mut_model &dm, mcont_t &, lookup_t &mut_lookup,
                 Args &&... args)
            -> decltype(std::bind(
                &discrete_mut_model::make_mut<traits::recycling_bin_t<mcont_t>,
                                              lookup_t, mcont_t>,
                &dm, std::placeholders::_1, std::placeholders::_2,
                std::forward<Args>(args)..., std::ref(mut_lookup)))
        {
            return std::bind(
                &discrete_mut_model::make_mut<traits::recycling_bin_t<mcont_t>,
                                              lookup_t, mcont_t>,
                &dm, std::placeholders::_1, std::placeholders::_2,
                std::forward<Args>(args)..., std::ref(mut_lookup));
        }

        /*! Return a vector of callables bount
         *  to KTfwd::extensions::discrete_mut_model::make_mut
         */
        template <typename mcont_t, typename lookup_t, class... Args>
        inline auto
        bind_vec_dmm(const std::vector<discrete_mut_model> &vdm,
                     mcont_t &mutations, lookup_t &mut_lookup,
                     const gsl_rng *r,
                     const std::vector<double> &neutral_mutrates,
                     const std::vector<double> &selected_mutrates,
                     Args &&... args)
            -> std::vector<decltype(
                bind_dmm(vdm[0], mutations, mut_lookup, r, neutral_mutrates[0],
                         selected_mutrates[0], std::forward<Args>(args)...))>
        {
            if (vdm.size() != neutral_mutrates.size()
                || vdm.size() != selected_mutrates.size())
                {
                    throw std::invalid_argument(
                        "container sizes must all be equal");
                }
            std::vector<decltype(
                bind_dmm(vdm[0], mutations, mut_lookup, r, neutral_mutrates[0],
                         selected_mutrates[0], std::forward<Args>(args)...))>
                rv;
            std::size_t i = 0;
            for (auto &&dm : vdm)
                {
                    rv.emplace_back(bind_dmm(
                        dm, mutations, mut_lookup, r, neutral_mutrates[i],
                        selected_mutrates[i], std::forward<Args>(args)...));
                    ++i;
                }
            return rv;
        }

        struct discrete_rec_model_data
        {
            std::vector<double> beg, end, weight;
            discrete_rec_model_data(const std::vector<double> &&b,
                                    const std::vector<double> &&e,
                                    const std::vector<double> &&w)
                : beg(std::move(b)), end(std::move(e)), weight(std::move(w))
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
            KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr lookup;
            void
            assign_weights()
            {
                if (data->weight.size())
                    lookup = KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(
                        gsl_ran_discrete_preproc(data->weight.size(),
                                                 data->weight.data()));
            }

          public:
            /*!
              \param __beg Region beginnings
              \param __end Region ends
              \param __weight Region weights
            */
            discrete_rec_model(const std::vector<double> __beg,
                               const std::vector<double> __end,
                               const std::vector<double> __weight)
                : data(new discrete_rec_model_data(
                      std::move(__beg), std::move(__end), std::move(__weight)))
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
            template <typename gamete_t, typename mcont_t>
            inline result_type
            operator()(const gsl_rng *r, const double recrate,
                       const gamete_t &, const gamete_t &,
                       const mcont_t &) const
            {
                assert(!(recrate == 0.
                         && (data->beg.empty() || data->end.empty())));
                auto nbreaks = gsl_ran_poisson(r, recrate);
                if (!nbreaks)
                    return {};

                result_type rv;
                for (unsigned i = 0; i < nbreaks; ++i)
                    {
                        size_t region = gsl_ran_discrete(r, lookup.get());
                        rv.push_back(gsl_ran_flat(r, data->beg[region],
                                                  data->end[region]));
                    }
                std::sort(rv.begin(), rv.end());
                rv.push_back(std::numeric_limits<double>::max());
                return rv;
            }
        };

        /*!
          Returns a function call bound to discrete_rec_model::operator().

          See unit test extensions_regionsTest.cc for example usage.
         */
        template <typename gcont_t, typename mcont_t>
        inline auto
        bind_drm(const discrete_rec_model &drm, const gcont_t &,
                 const mcont_t &, const gsl_rng *r, const double recrate)
            -> decltype(std::bind(&discrete_rec_model::operator() <
                                      typename gcont_t::value_type,
                                  mcont_t >, &drm, r, recrate,
                                  std::placeholders::_1, std::placeholders::_2,
                                  std::placeholders::_3))
        {
            return std::bind(
                &discrete_rec_model::operator() < typename gcont_t::value_type,
                mcont_t >, &drm, r, recrate, std::placeholders::_1,
                std::placeholders::_2, std::placeholders::_3);
        }

        /*! Returns a vector of function calls bound to
         *  discrete_rec_model::operator()
         */
        template <typename gcont_t, typename mcont_t>
        inline auto
        bind_vec_drm(const std::vector<discrete_rec_model> &vdrm,
                     const gcont_t &gametes, const mcont_t &mutations,
                     const gsl_rng *r, const std::vector<double> &recrates)
            -> std::vector<decltype(bind_drm(vdrm[0], gametes, mutations, r,
                                             recrates[0]))>
        {
            if (vdrm.size() != recrates.size())
                {
                    throw std::invalid_argument("unequal container sizes");
                }
            std::vector<decltype(
                bind_drm(vdrm[0], gametes, mutations, r, recrates[0]))>
                rv;
            static_assert(
                traits::is_rec_model<typename decltype(rv)::value_type,
                                     typename gcont_t::value_type,
                                     mcont_t>::value,
                "bound object must be a valid recombination model");
            std::size_t i = 0;
            for (auto &&drm : vdrm)
                {
                    rv.emplace_back(
                        bind_drm(drm, gametes, mutations, r, recrates[i++]));
                }
            return rv;
        }
    }
}
#endif
