#ifndef __FWDPP_SUGAR_MUTATION_INFSITES_HPP__
#define __FWDPP_SUGAR_MUTATION_INFSITES_HPP__

#include <fwdpp/internal/recycling.hpp>
#include <fwdpp/sugar/popgenmut.hpp>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <type_traits>
namespace KTfwd
{

    /*!
      \brief Generic function object implementing
      the infinitely-many sites mutation model for
      "standard" population-genetic scenarios
      \ingroup sugar
    */
    struct infsites
    {
        /*!
          Helper function.

          Generate a new mutation position and update the lookup table w/that
          position.
         */
        template <typename position_t, typename lookup_table_t>
        inline typename std::result_of<position_t()>::type
        generate_mut_pos(const position_t &posmaker,
                         lookup_table_t &lookup) const
        {
            static_assert(
                std::is_convertible<
                    typename std::result_of<position_t()>::type,
                    double>::value,
                "The return type of posmaker must be convertible to double.");
            auto pos = posmaker();
            while (lookup.find(pos) != lookup.end())
                {
                    pos = posmaker();
                }
            lookup.insert(pos);
            return pos;
        }

        /*!
          \param r A const gsl_rng *
          \param lookup A lookup table of mutation positions, see @ref
          md_md_policies for details.
          \param generation Generation when this mutation is happening
          \param neutral_mutation_rate Either the rate at which neutral
          variants arise (per gamete per generation), or something directly
          proportional to it
          \param selected_mutation_rate Either the rate at which non-neutral
          variants arise (per gamete per generation), or something directly
          proportional to it
          \param posmaker A policy that returns the position of the new
          mutation
          \param smaker A policy generating the selection coefficient/effect
          size associated with non-neutral variants
          \param hmaker A policy generating the dominance associated with
          non-neutral variants

          \note A mutation will be "selected" with probability
          selected_mutation_rate/(selected_mutation_rate +
          neutral_mutation_rate)
        */
        template <typename queue_t, typename mcont_t, typename lookup_table_t,
                  typename position_t, typename sdist_t, typename hdist_t>
        inline
            typename std::enable_if<std::is_same<typename mcont_t::value_type,
                                                 popgenmut>::value,
                                    std::size_t>::type
            operator()(queue_t &recycling_bin, mcont_t &mutations,
                       const gsl_rng *r, lookup_table_t &lookup,
                       const uint_t &generation,
                       const double &neutral_mutation_rate,
                       const double &selected_mutation_rate,
                       const position_t &posmaker, const sdist_t &smaker,
                       const hdist_t &hmaker) const
        {
            static_assert(std::is_same<typename mcont_t::value_type,
                                       KTfwd::popgenmut>::value,
                          "mcont_t::value_type must be KTfwd::popgenmut");
            // Establish position of new mutation
            auto pos = this->generate_mut_pos(posmaker, lookup);
            bool selected
                = (gsl_rng_uniform(r)
                   < selected_mutation_rate
                         / (neutral_mutation_rate + selected_mutation_rate));
            return fwdpp_internal::recycle_mutation_helper(
                recycling_bin, mutations, pos, (selected) ? smaker() : 0.,
                (selected) ? hmaker() : 0., generation);
        }

        /*!
          \brief Overload for different position distributions for neutral and
          non-neutral variants

          \param r A const gsl_rng *
          \param lookup A lookup table of mutation positions, see @ref
          md_md_policies for details.
          \param generation Generation when this mutation is happening
          \param neutral_mutation_rate Either the rate at which neutral
          variants arise (per gamete per generation), or something directly
          proportional to it
          \param selected_mutation_rate Either the rate at which non-neutral
          variants arise (per gamete per generation), or something directly
          proportional to it
          \param nposmaker A policy that returns the position of new neutral
          mutation
          \param sposmaker A policy that returns the position of new
          non-neutral mutation
          \param smaker A policy generating the selection coefficient/effect
          size associated with non-neutral variants
          \param hmaker A policy generating the dominance associated with
          non-neutral variants

          \note A mutation will be "selected" with probability
          selected_mutation_rate/(selected_mutation_rate +
          neutral_mutation_rate)
        */
        template <typename queue_t, typename mcont_t, typename lookup_table_t,
                  typename nposition_t, typename sposition_t, typename sdist_t,
                  typename hdist_t>
        inline
            typename std::enable_if<std::is_same<typename mcont_t::value_type,
                                                 popgenmut>::value,
                                    std::size_t>::type
            operator()(queue_t &recycling_bin, mcont_t &mutations,
                       const gsl_rng *r, lookup_table_t &lookup,
                       const uint_t &generation,
                       const double &neutral_mutation_rate,
                       const double &selected_mutation_rate,
                       const nposition_t &nposmaker,
                       const sposition_t &sposmaker, const sdist_t &smaker,
                       const hdist_t &hmaker) const
        {
            bool selected
                = gsl_rng_uniform(r)
                  < selected_mutation_rate
                        / (neutral_mutation_rate + selected_mutation_rate);
            if (selected)
                {
                    auto pos = this->generate_mut_pos(sposmaker, lookup);
                    return fwdpp_internal::recycle_mutation_helper(
                        recycling_bin, mutations, pos, smaker(), hmaker(),
                        generation);
                }
            // Establish position of new mutation
            auto pos = this->generate_mut_pos(nposmaker, lookup);
            return fwdpp_internal::recycle_mutation_helper(
                recycling_bin, mutations, pos, 0., 0., generation);
        }

        /*!
          \brief Overload for different position distributions for neutral and
          non-neutral variants

          \param r A const gsl_rng *
          \param lookup A lookup table of mutation positions, see @ref
          md_md_policies for details.
          \param generation Generation when this mutation is happening
          \param neutral_mutation_rate Either the rate at which neutral
          variants arise (per gamete per generation), or something directly
          proportional to it
          \param selected_mutation_rate Either the rate at which non-neutral
          variants arise (per gamete per generation), or something directly
          proportional to it
          \param nposmaker A policy that returns the position of new neutral
          mutation
          \param sposmaker A policy that returns the position of new
          non-neutral mutation
          \param smaker A policy generating the selection coefficient/effect
          size associated with non-neutral variants
          \param hmaker A policy generating the dominance associated with
          non-neutral variants

          \note A mutation will be "selected" with probability
          selected_mutation_rate/(selected_mutation_rate +
          neutral_mutation_rate)
        */
        template <typename queue_t, typename mcont_t, typename lookup_table_t,
                  typename nposition_t, typename sposition_t, typename sdist_t,
                  typename hdist_t>
        inline
            typename std::enable_if<std::is_same<typename mcont_t::value_type,
                                                 popgenmut>::value,
                                    std::size_t>::type
            operator()(queue_t &recycling_bin, mcont_t &mutations,
                       const gsl_rng *r, lookup_table_t &lookup,
                       const uint_t *generation,
                       const double &neutral_mutation_rate,
                       const double &selected_mutation_rate,
                       const nposition_t &nposmaker,
                       const sposition_t &sposmaker, const sdist_t &smaker,
                       const hdist_t &hmaker) const
        {
            return this->operator()(recycling_bin, mutations, r, lookup,
                                    *generation, neutral_mutation_rate,
                                    selected_mutation_rate, nposmaker,
                                    sposmaker, smaker, hmaker);
        }
        /*!
          \param r A const gsl_rng *
          \param lookup A lookup table of mutation positions, see @ref
          md_md_policies for details.
          \param generation Generation when this mutation is happening
          \param neutral_mutation_rate Either the rate at which neutral
          variants arise (per gamete per generation), or something directly
          proportional to it
          \param selected_mutation_rate Either the rate at which non-neutral
          variants arise (per gamete per generation), or something directly
          proportional to it
          \param posmaker A policy that returns the position of the new
          mutation
          \param smaker A policy generating the selection coefficient/effect
          size associated with non-neutral variants
          \param hmaker A policy generating the dominance associated with
          non-neutral variants

          \note A mutation will be "selected" with probability
          selected_mutation_rate/(selected_mutation_rate +
          neutral_mutation_rate)
        */
        template <typename queue_t, typename mcont_t, typename lookup_table_t,
                  typename position_t, typename sdist_t, typename hdist_t>
        inline
            typename std::enable_if<std::is_same<typename mcont_t::value_type,
                                                 popgenmut>::value,
                                    std::size_t>::type
            operator()(queue_t &recycling_bin, mcont_t &mutations,
                       const gsl_rng *r, lookup_table_t &lookup,
                       const uint_t *generation,
                       const double &neutral_mutation_rate,
                       const double &selected_mutation_rate,
                       const position_t &posmaker, const sdist_t &smaker,
                       const hdist_t &hmaker) const
        {
            return this->operator()(recycling_bin, mutations, r, lookup,
                                    *generation, neutral_mutation_rate,
                                    selected_mutation_rate, posmaker, smaker,
                                    hmaker);
        }

        /*!
          \param r A const gsl_rng *
          \param lookup A lookup table of mutation positions, see @ref
          md_md_policies for details.
          \param neutral_mutation_rate Either the rate at which neutral
          variants arise (per gamete per generation), or something directly
          proportional to it
          \param selected_mutation_rate Either the rate at which non-neutral
          variants arise (per gamete per generation), or something directly
          proportional to it
          \param posmaker A policy that returns the position of the new
          mutation
          \param smaker A policy generating the selection coefficient/effect
          size associated with non-neutral variants
          \param hmaker A policy generating the dominance associated with
          non-neutral variants

          \note A mutation will be "selected" with probability
          selected_mutation_rate/(selected_mutation_rate +
          neutral_mutation_rate)
        */
        template <typename queue_t, typename mcont_t, typename lookup_table_t,
                  typename position_t, typename sdist_t, typename hdist_t>
        inline
            typename std::enable_if<std::is_same<typename mcont_t::value_type,
                                                 mutation>::value,
                                    std::size_t>::type
            operator()(queue_t &mutation_recycling_bin, mcont_t &mutations,
                       const gsl_rng *r, lookup_table_t &lookup,
                       const double &neutral_mutation_rate,
                       const double &selected_mutation_rate,
                       const position_t &posmaker, const sdist_t &smaker,
                       const hdist_t &hmaker) const
        {
            static_assert(std::is_same<typename mcont_t::value_type,
                                       KTfwd::mutation>::value,
                          "mcont_t::value_type must be KTfwd::mutation");
            // Establish position of new mutation
            auto pos = this->generate_mut_pos(posmaker, lookup);
            bool selected
                = (gsl_rng_uniform(r)
                   < selected_mutation_rate
                         / (neutral_mutation_rate + selected_mutation_rate));
            return fwdpp_internal::recycle_mutation_helper(
                mutation_recycling_bin, mutations, pos,
                (selected) ? smaker() : 0., (selected) ? hmaker() : 0.);
        }

        /*!
          \brief Overload for different position distributions for neutral and
          non-neutral variants

          \param r A const gsl_rng *
          \param lookup A lookup table of mutation positions, see @ref
          md_md_policies for details.
          \param neutral_mutation_rate Either the rate at which neutral
          variants arise (per gamete per generation), or something directly
          proportional to it
          \param selected_mutation_rate Either the rate at which non-neutral
          variants arise (per gamete per generation), or something directly
          proportional to it
          \param nposmaker A policy that returns the position of new neutral
          mutation
          \param sposmaker A policy that returns the position of new
          non-neutral mutation
          \param smaker A policy generating the selection coefficient/effect
          size associated with non-neutral variants
          \param hmaker A policy generating the dominance associated with
          non-neutral variants

          \note A mutation will be "selected" with probability
          selected_mutation_rate/(selected_mutation_rate +
          neutral_mutation_rate)
        */
        template <typename queue_t, typename mcont_t, typename lookup_table_t,
                  typename nposition_t, typename sposition_t, typename sdist_t,
                  typename hdist_t>
        inline
            typename std::enable_if<std::is_same<typename mcont_t::value_type,
                                                 mutation>::value,
                                    std::size_t>::type
            operator()(queue_t &recycling_bin, mcont_t &mutations,
                       const gsl_rng *r, lookup_table_t &lookup,
                       const double &neutral_mutation_rate,
                       const double &selected_mutation_rate,
                       const nposition_t &nposmaker,
                       const sposition_t &sposmaker, const sdist_t &smaker,
                       const hdist_t &hmaker) const
        {
            bool selected
                = gsl_rng_uniform(r)
                  < selected_mutation_rate
                        / (neutral_mutation_rate + selected_mutation_rate);
            if (selected)
                {
                    auto pos = this->generate_mut_pos(sposmaker, lookup);
                    return fwdpp_internal::recycle_mutation_helper(
                        recycling_bin, mutations, pos, smaker(), hmaker());
                }
            auto pos = this->generate_mut_pos(nposmaker, lookup);
            // return a neutral mutation
            return fwdpp_internal::recycle_mutation_helper(
                recycling_bin, mutations, pos, 0., 1, 0.);
        }
    };
}

#endif
