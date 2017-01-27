#ifndef __FWDPP_SUGAR_MULTILOC_MULTILOC_HPP__
#define __FWDPP_SUGAR_MULTILOC_MULTILOC_HPP__

#include <type_traits>
#include <vector>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/sugar/poptypes/tags.hpp>
#include <fwdpp/sugar/poptypes/popbase.hpp>

namespace KTfwd
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
        struct multiloc
            : public popbase<mutation_type, mcont, gcont, dipvector, mvector,
                             ftvector, lookup_table_type>
        {
            //! Dispatch tags for other parts of sugar layer
            using popmodel_t = sugar::MULTILOCPOP_TAG;
            //! Typedef for base class
            using popbase_t = popbase<mutation_type, mcont, gcont, dipvector,
                                      mvector, ftvector, lookup_table_type>;

            //! Population size
            uint_t N;

            //! Container of individuals
            typename popbase_t::dipvector_t diploids;

            //! Construct with population size and number of loci
            multiloc(
                const uint_t &__N, const uint_t &__nloci,
                typename popbase_t::gamete_t::mutation_container::size_type
                    reserve_size
                = 100)
                : popbase_t(__nloci * __N, reserve_size), N(__N),
                  diploids(__N, typename popbase_t::diploid_t(
                                    __nloci,
                                    typename popbase_t::diploid_t::value_type(
                                        0, 0)))
            {
            }

            bool
            operator==(const multiloc &rhs) const
            {
                return this->diploids == rhs.diploids
                       && popbase_t::is_equal(rhs);
            }

            //! Empty all containers
            void
            clear()
            {
                diploids.clear();
                popbase_t::clear_containers();
            }
        };
    }
}

#endif
