#ifndef __FWDPP_SUGAR_METAPOP_METAPOP_HPP__
#define __FWDPP_SUGAR_METAPOP_METAPOP_HPP__

#include <type_traits>
#include <numeric>
#include <fwdpp/sugar/poptypes/singlepop.hpp>
#include <fwdpp/sugar/poptypes/tags.hpp>
#include <fwdpp/sugar/poptypes/popbase.hpp>

namespace KTfwd
{
    namespace sugar
    {
        /*!
          \brief Abstraction of what is needed to simulate a metapopulation
          using an individual-based sampler from fwdpp

          All that is missing is the mutation_type and the container types.

          See @ref md_md_sugar for rationale, etc.

          \ingroup sugar
        */
        template <typename mutation_type, typename mcont, typename gcont,
                  typename dipvector, typename vdipvector, typename mvector,
                  typename ftvector, typename lookup_table_type>
        class metapop : public popbase<mutation_type, mcont, gcont, dipvector,
                                       mvector, ftvector, lookup_table_type>
        {
          private:
            void
            init_vectors()
            {
                for (uint_t i = 0; i < Ns.size(); ++i)
                    {
                        diploids.emplace_back(typename popbase_t::dipvector_t(
                            Ns[i], typename dipvector::value_type(0, 0)));
                    }
            }

          public:
            virtual ~metapop() = default;
			metapop(metapop &&) = default;
			metapop(const metapop &) = default;
            metapop &operator=(metapop &&) = default;
            metapop &operator=(const metapop &) = default;
            //! Deme sizes
            std::vector<uint_t> Ns;

            //! Typedef for base class
            using popbase_t = popbase<mutation_type, mcont, gcont, dipvector,
                                      mvector, ftvector, lookup_table_type>;

            //! Dispatch tags for other parts of sugar layer
            using popmodel_t = sugar::METAPOP_TAG;

            //! Typedef for vector<vector<popbase_t::diploid_t>
            using vdipvector_t = vdipvector;

            //! Fitness function signature compatible with this type
            using fitness_t
                = KTfwd::traits::fitness_fxn_t<typename popbase_t::dipvector_t,
                                               typename popbase_t::gcont_t,
                                               typename popbase_t::mcont_t>;
            //! Metapops can be constructed from singlepops of this type
            using compat_singlepop_t = sugar::singlepop<
                typename popbase_t::mutation_t, typename popbase_t::mcont_t,
                typename popbase_t::gcont_t, typename popbase_t::dipvector_t,
                typename popbase_t::mvector_t, typename popbase_t::ftvector_t,
                typename popbase_t::lookup_table_t>;
            //! Container of diploids
            vdipvector_t diploids;

            //! Construct with a cont of deme sizes
            explicit metapop(
                std::initializer_list<uint_t> __Ns,
                typename popbase_t::gamete_t::mutation_container::size_type
                    reserve_size
                = 100)
                : popbase_t(
                      std::accumulate(std::begin(__Ns), std::end(__Ns), 0.),
                      reserve_size),
                  Ns(__Ns), diploids(vdipvector_t())
            {
                init_vectors();
            }

            //! Construct with array of deme sizes
            explicit metapop(
                const uint_t *__Ns, const size_t num_Ns,
                typename popbase_t::gamete_t::mutation_container::size_type
                    reserve_size
                = 100)
                : popbase_t(std::accumulate(__Ns, __Ns + num_Ns, 0.),
                            reserve_size),
                  Ns(std::vector<uint_t>(__Ns, __Ns + num_Ns)),
                  diploids(vdipvector_t())
            {
                init_vectors();
            }

            //! Copy construct from a singlepop based on the same basic types
            explicit metapop(const compat_singlepop_t &spop)
                : popbase_t(0), Ns({ spop.N })
            {
                this->mutations = spop.mutations;
                this->mcounts = spop.mcounts;
                this->gametes = spop.gametes;
                this->diploids = vdipvector_t(1, spop.diploids);
                this->neutral = spop.neutral;
                this->selected = spop.selected;
                this->mut_lookup = spop.mut_lookup;
                this->fixations = spop.fixations;
                this->fixation_times = spop.fixation_times;
            }

            //! Move construct from a singlepop based on the same basic types
            metapop(compat_singlepop_t &&spop) : popbase_t(0), Ns({ spop.N })
            {
                this->mutations = std::move(spop.mutations);
                this->mcounts = std::move(spop.mcounts);
                this->gametes = std::move(spop.gametes);
                this->diploids = vdipvector_t(1, std::move(spop.diploids));
                this->neutral = std::move(spop.neutral);
                this->selected = std::move(spop.selected);
                this->mut_lookup = std::move(spop.mut_lookup);
                this->fixations = std::move(spop.fixations);
                this->fixation_times = std::move(spop.fixation_times);
            }

            bool
            operator==(const metapop &rhs) const
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
