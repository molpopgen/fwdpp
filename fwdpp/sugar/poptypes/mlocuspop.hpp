#ifndef __FWDPP_SUGAR_MULTILOC_MULTILOC_HPP__
#define __FWDPP_SUGAR_MULTILOC_MULTILOC_HPP__

#include <type_traits>
#include <vector>
#include <utility>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/sugar/poptypes/tags.hpp>
#include <fwdpp/sugar/poptypes/popbase.hpp>

namespace fwdpp
{
    namespace sugar
    {
        /*!
          \brief Abstraction of what is needed to simulate a multilocus
          simulation using an individual-based sampler from fwdpp.

          All that is missing is the mutation_type and the container types.

          See @ref md_md_sugar for rationale, etc.

          \ingroup sugar
        */
        template <typename mutation_type, typename mcont, typename gcont,
                  typename dipvector, typename mvector, typename ftvector,
                  typename lookup_table_type>
        class mlocuspop : public popbase<mutation_type, mcont, gcont, mvector,
                                         ftvector, lookup_table_type>
        {
          private:
            void
            process_individual_input()
            {
                std::vector<uint_t> gcounts(this->gametes.size(), 0);
                for (auto &&locus : diploids)
                    {
                        for (auto &&dip : locus)
                            {
                                this->validate_individual_keys(dip.first);
                                this->validate_individual_keys(dip.second);
                                gcounts[dip.first]++;
                                gcounts[dip.second]++;
                            }
                    }
                this->validate_gamete_counts(gcounts);
            }

          public:
            virtual ~mlocuspop() = default;
            mlocuspop(mlocuspop &&) = default;
            mlocuspop(const mlocuspop &) = default;
            mlocuspop &operator=(mlocuspop &&) = default;
            mlocuspop &operator=(const mlocuspop &) = default;
            using dipvector_t = dipvector;
            using diploid_t = typename dipvector::value_type;
            //! Dispatch tags for other parts of sugar layer
            using popmodel_t = sugar::MULTILOC_TAG;
            //! Typedef for base class
            using popbase_t = popbase<mutation_type, mcont, gcont, mvector,
                                      ftvector, lookup_table_type>;

            static_assert(traits::is_multilocus_diploid_t<diploid_t>::value,
                          "Require that dipvector_t::value_type be a diploid");

            //! Population size
            uint_t N;

            //! Container of individuals
            dipvector_t diploids;

            /*! The positional boundaries of each locus/region,
             *  expressed as half-open intervals [min,max).
             */
            std::vector<std::pair<double, double>> locus_boundaries;

            /*! Construct with population size, number of loci,
             *  and locus boundaries.
             */
            mlocuspop(
                const uint_t &__N, const uint_t &__nloci,
                const std::vector<std::pair<double, double>> &locus_boundaries_
                = std::vector<std::pair<double, double>>(),
                typename popbase_t::gamete_t::mutation_container::size_type
                    reserve_size
                = 100)
                : popbase_t(__nloci * 2 * __N, reserve_size), N(__N),
                  diploids(__N, diploid_t(__nloci, { 0, 0 })),
                  locus_boundaries(locus_boundaries_)
            {
            }

            template <typename diploids_input, typename gametes_input,
                      typename mutations_input>
            explicit mlocuspop(
                diploids_input &&d, gametes_input &&g, mutations_input &&m,
                const std::vector<std::pair<double, double>> &locus_boundaries_
                = std::vector<std::pair<double, double>>())
                : popbase_t(std::forward<gametes_input>(g),
                            std::forward<mutations_input>(m), 100),
                  N(d.size()), diploids(std::forward<diploids_input>(d)),
                  locus_boundaries(locus_boundaries_)
            {
                this->process_individual_input();
            }

            bool
            operator==(const mlocuspop &rhs) const
            {
                return this->diploids == rhs.diploids
                       && this->locus_boundaries == rhs.locus_boundaries
                       && popbase_t::is_equal(rhs);
            }

            //! Empty all containers
            void
            clear()
            {
                diploids.clear();
                locus_boundaries.clear();
                popbase_t::clear_containers();
            }
        };
    } // namespace sugar
} // namespace fwdpp

#endif
