#ifndef FWDPP_SUGAR_POPBASE_HPP
#define FWDPP_SUGAR_POPBASE_HPP

#include <type_traits>
#include <vector>
#include <exception>
#include <fwdpp/type_traits.hpp>
#include <fwdpp/internal/sample_diploid_helpers.hpp>

namespace KTfwd
{
    namespace sugar
    {
        template <typename mutation_type, typename mcont, typename gcont,
                  typename dipvector, typename mvector, typename ftvector,
                  typename lookup_table_type>
        class popbase
        /*!
          \ingroup sugar
          \brief Base class for population objects
          \note Added in fwdpp 0.5.0
         */
        {
            static_assert(typename KTfwd::traits::is_gamete<
                              typename gcont::value_type>::type(),
                          "gcont::value_type must be a gamete type");
            static_assert(typename KTfwd::traits::is_mutation<
                              typename mcont::value_type>::type(),
                          "mcont::value_type must be a mutation type");
            /// This function is used to validate input
            /// data.  It should be called from constructors
            /// allowing the initialization of populations
            /// from pre-calculated diploids, gametes, and mutations.
            /// \version
            /// Added in 0.5.7
            virtual void process_diploid_input() = 0;
            void check_mutation_keys(
                const typename gcont::value_type::mutation_container &m,
                const mcont &mutations, const bool neutrality);
            void fill_internal_structures();
          protected:
            // Protected members introduced in 0.5.7 to help
            // derived classes implement constructors based
            // on user input of population data.
            void
            validate_diploid_keys(const std::size_t first,
                                  const std::size_t second)
            {
                if (first >= this->gametes.size()
                    || second >= this->gametes.size())
                    {
                        throw std::out_of_range(
                            "diploid contains out of range keys");
                    }
                if (!this->gametes[first].n || !this->gametes[second].n)
                    {
                        throw std::runtime_error("diploid refers to "
                                                 "gamete marked as "
                                                 "extinct");
                    }
            }
            void
            validate_gamete_counts(const std::vector<uint_t> &gcounts)
            {
                for (std::size_t i = 0; i < gcounts.size(); ++i)
                    {
                        if (gcounts[i] != this->gametes[i].n)
                            {
                                throw std::runtime_error(
                                    "gamete count does not match number of "
                                    "diploids referring to it");
                            }
                    }
            }

          public:
            virtual ~popbase() = default;
            popbase(popbase &&) = default;
            popbase(const popbase &) = default;
            popbase &operator=(popbase &&) = default;
            popbase &operator=(const popbase &) = default;
            //! Mutation type
            using mutation_t = mutation_type;
            //! Gamete type
            using gamete_t = typename gcont::value_type;
            //! Diploid vector type
            using dipvector_t = dipvector;
            //! Diploid type
            using diploid_t = typename dipvector_t::value_type;
            //! Mutation vec type
            using mcont_t = mcont;
            //! Mutation count vector type
            using mcount_t = std::vector<uint_t>;
            //! Gamete vec type
            using gcont_t = gcont;
            //! Lookup table type for recording mutation positions, etc.
            using lookup_table_t = lookup_table_type;
            //! container type for fixations
            using mvector_t = mvector;
            //! container type for fixation times
            using ftvector_t = ftvector;

            mcont_t mutations;
            /*!
              Used to keep track of mutation frequencies.

              Should have memory reserved externally,
              based on some good guess.
            */

            //! Contains number of times each mutation exists
            mcount_t mcounts;
            //! Container of gametes
            gcont gametes;

            /*!
                          Containers that can be used as intermediates during
              the generation
                          of new gametes.

              The requirement to declare these was introduced in fwdpp 0.3.3.

              In previous versions of the library, vectors like this had to be
              allocated
              for every crossover event for every generation.  The result was
              an excessive
              number of requests for memory allocation.

              Now, we create the vector once per replicate.  Further, we will
              reserve memory
              here, to minimize reallocs, etc., within fwdpp.

              Internally, fwdpp's job is to make sure that this vector is
              appropriately
              and efficiently cleared, but only when needed.

              \note: if not using the sugar features, you can create these
              vectors
              only once per simulation...
            */
            typename gamete_t::mutation_container neutral, selected;

            /*!
              \brief Can be used to track positions of segregating mutations.
              \note Must have interface like std::map or std::unordered_set
            */
            lookup_table_type mut_lookup;
            //! Vector of mutation_t to track fixations
            mvector fixations;
            /*! \brief vector<uint_t> records times when mutation_ts
              were added to mut_lookup
            */
            ftvector fixation_times;

            //! Constructor
            explicit popbase(
                const uint_t &popsize,
                typename gamete_t::mutation_container::size_type reserve_size
                = 100)
                : // No muts in the population
                  mutations(mcont_t()),
                  mcounts(mcount_t()),
                  // The population contains a single gamete in 2N copies
                  gametes(gcont(1, gamete_t(2 * popsize))),
                  neutral(typename gamete_t::mutation_container()),
                  selected(typename gamete_t::mutation_container()),
                  mut_lookup(lookup_table_type()), fixations(mvector()),
                  fixation_times(ftvector())
            {
                // This is a good number for reserving,
                // allowing for extra allocations when recycling is doing its
                // thing
                gametes.reserve(4 * popsize);
                // Reserve memory
                neutral.reserve(reserve_size);
                selected.reserve(reserve_size);
            }

            template <typename gametes_input, typename mutations_input>
            explicit popbase(
                gametes_input &&g, mutations_input &m,
                typename gamete_t::mutation_container::size_type reserve_size)
                : mutations(std::forward<mutations_input>(m)), mcounts{},
                  gametes(std::forward<gametes_input>(g)), neutral{},
                  selected{}, mut_lookup{}, fixations{}, fixation_times{}
            {
                this->fill_internal_structures();
                neutral.reserve(reserve_size);
                selected.reserve(reserve_size);
            }

            bool
            is_equal(const popbase &rhs) const
            {
                return this->mutations == rhs.mutations
                       && this->mcounts == rhs.mcounts
                       && this->gametes == rhs.gametes
                       && this->fixations == rhs.fixations
                       && this->fixation_times == rhs.fixation_times;
            }

            //! Empty all the containers
            void
            clear_containers()
            {
                mutations.clear();
                mcounts.clear();
                gametes.clear();
                mut_lookup.clear();
                fixations.clear();
                fixation_times.clear();
            }
        };

        template <typename mutation_type, typename mcont, typename gcont,
                  typename dipvector, typename mvector, typename ftvector,
                  typename lookup_table_type>
        void
        popbase<mutation_type, mcont, gcont, dipvector, mvector, ftvector,
                lookup_table_type>::
            check_mutation_keys(
                const typename gcont::value_type::mutation_container &m,
                const mcont &mutations, const bool neutrality)
        {
            for (const auto &k : m)
                {
                    mcounts.resize(mutations.size(), 0);
                    if (k >= mutations.size())
                        {
                            throw std::out_of_range(
                                "gamete contains mutation key that is "
                                "out of range");
                        }
                    if (mutations[k].neutral != neutrality)
                        {
                            throw std::logic_error(
                                "gamete contains key to "
                                "mutation in wrong container.");
                        }
                }
        }

        template <typename mutation_type, typename mcont, typename gcont,
                  typename dipvector, typename mvector, typename ftvector,
                  typename lookup_table_type>
        void
        popbase<mutation_type, mcont, gcont, dipvector, mvector, ftvector,
                lookup_table_type>::fill_internal_structures()
        {
            mut_lookup.clear();
            fixations.clear();
            fixation_times.clear();
            mcounts.clear();
            for (auto &&m : mutations)
                {
                    mut_lookup.insert(m.pos);
                }
            for (const auto &g : gametes)
                {
                    check_mutation_keys(g.mutations, mutations, true);
                    check_mutation_keys(g.smutations, mutations, false);
                }
            fwdpp_internal::process_gametes(gametes, mutations, mcounts);
        }
    }
}
#endif
