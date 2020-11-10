#ifndef __FWDPP_SUGAR_SINGLEPOP_SINGLEPOP_HPP__
#define __FWDPP_SUGAR_SINGLEPOP_SINGLEPOP_HPP__

/*
  A structure representing a single-locus population.
  The user initizializes it with a population size, N
*/

#include <fwdpp/poptypes/tags.hpp>
#include <fwdpp/poptypes/popbase.hpp>

namespace fwdpp
{
    namespace poptypes
    {
        /*!
          \brief Abstraction of what is needed to simulate a
          single-locus population.

          All that is missing is the mutation_type and the container types.
        */
        template <typename mutation_type, typename mcont, typename gcont,
                  typename dipvector, typename mvector, typename ftvector,
                  typename lookup_table_type>
        class diploid_population : public popbase<mutation_type, mcont, gcont, mvector,
                                                  ftvector, lookup_table_type>
        {
          public:
            virtual ~diploid_population() = default;
            diploid_population(diploid_population &&) = default;
            diploid_population(const diploid_population &) = default;
            diploid_population &operator=(diploid_population &&) = default;
            diploid_population &operator=(const diploid_population &) = default;
            //! Population size
            uint_t N;

            using dipvector_t = dipvector;
            using diploid_type = typename dipvector::value_type;
            using diploid_t [[deprecated("use diploid_type")]] =
                typename dipvector::value_type;
            //! Typedef for base class
            using popbase_t = popbase<mutation_type, mcont, gcont, mvector, ftvector,
                                      lookup_table_type>;
            //! Dispatch tag
            using popmodel_t = poptypes::DIPLOID_TAG;
            //! Fitness function signature compatible with this type
            using fitness_t
                = fwdpp::traits::fitness_fxn_t<dipvector_t,
                                               typename popbase_t::genome_container,
                                               typename popbase_t::mutation_container>;

            //! Container of diploids
            dipvector_t diploids;

            //! Constructor
            explicit diploid_population(
                const uint_t &popsize,
                typename popbase_t::haploid_genome_t::mutation_container::size_type
                    reserve_size
                = 100)
                : popbase_t(2 * popsize, reserve_size), N(popsize),
                  // All N diploids contain the only haploid_genome in the pop
                  diploids(dipvector_t(popsize, diploid_t(0, 0)))
            {
            }

            bool
            operator==(const diploid_population &rhs) const
            {
                return this->diploids == rhs.diploids && popbase_t::is_equal(rhs);
            }

            //! Empty all the containers
            void
            clear()
            {
                diploids.clear();
                popbase_t::clear_containers();
            }
        };
    } // namespace poptypes
} // namespace fwdpp
#endif
