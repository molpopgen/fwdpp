#ifndef __FWDPP_SUGAR_SINGLEPOP_SINGLEPOP_HPP__
#define __FWDPP_SUGAR_SINGLEPOP_SINGLEPOP_HPP__

/*
  A structure representing a single Wright-Fisher population.
  The user initizializes it with a population size, N
*/

#include <fwdpp/sugar/poptypes/tags.hpp>
#include <fwdpp/sugar/poptypes/popbase.hpp>
namespace KTfwd
{
    namespace sugar
    {
        /*!
          \brief Abstraction of what is needed to simulate a single population
          using an individual-based sampler from fwdpp

          All that is missing is the mutation_type and the container types.

          See @ref md_md_sugar for rationale, etc.

          \ingroup sugar
        */
        template <typename mutation_type, typename mcont, typename gcont,
                  typename dipvector, typename mvector, typename ftvector,
                  typename lookup_table_type>
        class singlepop
            : public popbase<mutation_type, mcont, gcont, dipvector, mvector,
                             ftvector, lookup_table_type>
        {
          public:
            virtual ~singlepop() = default;
            singlepop(singlepop &&) = default;
            singlepop(const singlepop &) = default;
            singlepop &operator=(singlepop &&) = default;
            singlepop &operator=(const singlepop &) = default;
            //! Population size
            uint_t N;

            //! Typedef for base class
            using popbase_t = popbase<mutation_type, mcont, gcont, dipvector,
                                      mvector, ftvector, lookup_table_type>;
            //! Dispatch tag for other parts of sugar layer
            using popmodel_t = sugar::SINGLEPOP_TAG;
            //! Fitness function signature compatible with this type
            using fitness_t
                = KTfwd::traits::fitness_fxn_t<typename popbase_t::dipvector_t,
                                               typename popbase_t::gcont_t,
                                               typename popbase_t::mcont_t>;

            //! Container of diploids
            typename popbase_t::dipvector_t diploids;

            //! Constructor
            explicit singlepop(
                const uint_t &popsize,
                typename popbase_t::gamete_t::mutation_container::size_type
                    reserve_size
                = 100)
                : popbase_t(popsize, reserve_size), N(popsize),
                  // All N diploids contain the only gamete in the pop
                  diploids(typename popbase_t::dipvector_t(
                      popsize, typename popbase_t::diploid_t(0, 0)))
            {
            }

            bool
            operator==(const singlepop &rhs) const
            {
                return this->diploids == rhs.diploids
                       && popbase_t::is_equal(rhs);
            }

            //! Empty all the containers
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
