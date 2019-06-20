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
        class diploid_population
            : public popbase<mutation_type, mcont, gcont, mvector, ftvector,
                             lookup_table_type>
        {
          private:
            void
            process_individual_input()
            {
                std::vector<uint_t> gcounts(this->haploid_genomes.size(), 0);
                for (auto &&dip : diploids)
                    {
                        this->validate_individual_keys(dip.first);
                        this->validate_individual_keys(dip.second);
                        gcounts[dip.first]++;
                        gcounts[dip.second]++;
                    }
                this->validate_haploid_genome_counts(gcounts);
            }

          public:
            virtual ~diploid_population() = default;
            diploid_population(diploid_population &&) = default;
            diploid_population(const diploid_population &) = default;
            diploid_population &operator=(diploid_population &&) = default;
            diploid_population &operator=(const diploid_population &)
                = default;
            //! Population size
            uint_t N;

            using dipvector_t = dipvector;
            using diploid_t = typename dipvector::value_type;
            //! Typedef for base class
            using popbase_t = popbase<mutation_type, mcont, gcont, mvector,
                                      ftvector, lookup_table_type>;
            //! Dispatch tag
            using popmodel_t = poptypes::DIPLOID_TAG;
            //! Fitness function signature compatible with this type
            using fitness_t
                = fwdpp::traits::fitness_fxn_t<dipvector_t,
                                               typename popbase_t::gcont_t,
                                               typename popbase_t::mcont_t>;

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

            template <typename diploids_input, typename haploid_genomes_input,
                      typename mutations_input>
            explicit diploid_population(diploids_input &&d, haploid_genomes_input &&g,
                                        mutations_input &&m)
                : popbase_t(std::forward<haploid_genomes_input>(g),
                            std::forward<mutations_input>(m), 100),
                  N{ static_cast<decltype(N)>(d.size()) },
                  diploids(std::forward<diploids_input>(d))
            //! Constructor for pre-determined population status
            {
                this->process_individual_input();
            }

            bool
            operator==(const diploid_population &rhs) const
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
    } // namespace poptypes
} // namespace fwdpp
#endif
