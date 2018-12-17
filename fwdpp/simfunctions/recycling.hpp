#ifndef FWDPP_RECYCLING
#define FWDPP_RECYCLING

#include <queue>
#include <stdexcept>
#include <type_traits>
#include <fwdpp/named_type.hpp>

namespace fwdpp
{
    namespace tags
    {
        struct mutation_recycling
        {
        };
        struct gamete_recycling
        {
        };
    } // namespace tags

    /// \brief FIFO queue for mutation recycling
    /// \version 0.7.4 added to fwdpp
    using flagged_mutation_queue
        = strong_types::named_type<std::queue<std::size_t>,
                                   tags::mutation_recycling>;

    /// \brief FIFO queue for gamete recycling
    /// \version 0.7.4 added to fwdpp
    using flagged_gamete_queue
        = strong_types::named_type<std::queue<std::size_t>,
                                   tags::gamete_recycling>;

    inline flagged_mutation_queue
    empty_mutation_queue()
    /// \brief Generate an empty flagged_mutation_queue
    {
        return flagged_mutation_queue(std::queue<std::size_t>());
    }

    inline flagged_gamete_queue
    empty_gamete_queue()
    //// \brief Generate an empty flagged_gamete_queue
    {
        return flagged_gamete_queue(std::queue<std::size_t>());
    }

    template <typename mcount_vec>
    inline flagged_mutation_queue
    make_mut_queue(const mcount_vec &mcounts)
    /// \brief Make a FIFO recycling queue for mutations
    /// \param mcounts Vector of mutation counts
    /// \note Simulations with tree sequences should use fwdpp::ts::make_mut_queue
    {
        std::queue<std::size_t> rv;
        const auto msize = mcounts.size();
        for (typename mcount_vec::size_type i = 0; i < msize; ++i)
            {
                if (!mcounts[i])
                    rv.push(i);
            }
        return flagged_mutation_queue(std::move(rv));
    }

    template <typename gvec_t>
    inline flagged_gamete_queue
    make_gamete_queue(const gvec_t &gametes)
    /// \brief Make a FIFO queue for recycling extinct gametes
    /// \param gametes Vector of gametes
    {
        std::queue<std::size_t> rv;
        const auto gsize = gametes.size();
        for (typename gvec_t::size_type i = 0; i < gsize; ++i)
            {
                if (!gametes[i].n)
                    rv.push(i);
            }
        return flagged_gamete_queue(std::move(rv));
    }

    template <typename gcont_t>
    inline std::size_t
    recycle_gamete(gcont_t &gametes,
                   flagged_gamete_queue &gamete_recycling_bin,
                   typename gcont_t::value_type::mutation_container &neutral,
                   typename gcont_t::value_type::mutation_container &selected)
    /// \brief Return location of a new gamete, recycling available memory if possible
    /// \param gametes vector of gametes
    /// \param gamete_recycling_bin A flagged_gamete_queue
    /// \param neutral Data for new gamete's neutral variants
    /// \param selected Data for new gamete's selected variants
    /// \return A location in \a gametes
    {
        // Try to recycle
        auto &ref = gamete_recycling_bin.get();
        if (!ref.empty())
            {
                auto idx = ref.front();
                ref.pop();
#ifndef NDEBUG
                if (gametes[idx].n)
                    {
                        throw std::runtime_error(
                            "FWDPP DEBUG: attempting to recycle an extant "
                            "gamete");
                    }
#endif
                gametes[idx].mutations.swap(neutral);
                gametes[idx].smutations.swap(selected);
                return idx;
            }
        gametes.emplace_back(0u, std::move(neutral), std::move(selected));
        return (gametes.size() - 1);
    }

    /*!
          \brief Helper function for mutation policies

          This function minimizes code duplication when writing mutation
          models.  It abstracts
          the operations needed to recycle an extinct mutation.

          \param mutation_recycling_bin  A FIFO queue of iterators pointing to
          extinct mutations.
          \param mutations A list of mutation objects
          \param args Parameter pack to be passed to constructor of an
          mcont_t::value_type
         */
    template <typename mcont_t, class... Args>
    inline std::size_t
    recycle_mutation_helper(flagged_mutation_queue &mutation_recycling_bin,
                            mcont_t &mutations, Args &&... args)
    /// \brief Helper function for implementing mutation generation functions.
    /// \param mutation_recycling_bin A flagged_mutation_queue
    /// \param mutations Container of mutations
    /// \param args Constructor arguments to create a new mutation
    /// \returns the location of the new variant in \a mutations
    {
        auto &ref = mutation_recycling_bin.get();
        if (!ref.empty())
            {
                auto rv = ref.front();
                ref.pop();
                mutations[rv] =
                    typename mcont_t::value_type(std::forward<Args>(args)...);
                return rv;
            }
        mutations.emplace_back(std::forward<Args>(args)...);
        return mutations.size() - 1;
    }
} // namespace fwdpp

#endif
