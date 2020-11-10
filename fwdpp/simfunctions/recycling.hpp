#ifndef FWDPP_RECYCLING
#define FWDPP_RECYCLING

#include <queue>
#include <stdexcept>
#include <type_traits>
#include <fwdpp/util/named_type.hpp>

namespace fwdpp
{
    namespace tags
    {
        struct mutation_recycling
        {
        };
        struct haploid_genome_recycling
        {
        };
    } // namespace tags

    /// \brief FIFO queue for mutation recycling
    /// \version 0.7.4 added to fwdpp
    using flagged_mutation_queue
        = strong_types::named_type<std::queue<std::size_t>, tags::mutation_recycling>;

    /// \brief FIFO queue for haploid_genome recycling
    /// \version 0.7.4 added to fwdpp
    using flagged_haploid_genome_queue
        = strong_types::named_type<std::queue<std::size_t>,
                                   tags::haploid_genome_recycling>;

    inline flagged_mutation_queue
    empty_mutation_queue()
    /// \brief Generate an empty flagged_mutation_queue
    {
        return flagged_mutation_queue(flagged_mutation_queue::value_type());
    }

    inline flagged_haploid_genome_queue
    empty_haploid_genome_queue()
    //// \brief Generate an empty flagged_haploid_genome_queue
    {
        return flagged_haploid_genome_queue(flagged_haploid_genome_queue::value_type());
    }

    template <typename mcount_vec>
    inline flagged_mutation_queue
    make_mut_queue(const mcount_vec &mcounts)
    /// \brief Make a FIFO recycling queue for mutations
    /// \param mcounts Vector of mutation counts
    /// \note Simulations with tree sequences should use fwdpp::ts::make_mut_queue
    {
        flagged_mutation_queue::value_type rv;
        const auto msize = mcounts.size();
        for (typename mcount_vec::size_type i = 0; i < msize; ++i)
            {
                if (!mcounts[i])
                    rv.push(i);
            }
        return flagged_mutation_queue(std::move(rv));
    }

    template <typename gvec_t>
    inline flagged_haploid_genome_queue
    make_haploid_genome_queue(const gvec_t &haploid_genomes)
    /// \brief Make a FIFO queue for recycling extinct haploid_genomes
    /// \param haploid_genomes Vector of haploid_genomes
    {
        flagged_haploid_genome_queue::value_type rv;
        const auto gsize = haploid_genomes.size();
        for (typename gvec_t::size_type i = 0; i < gsize; ++i)
            {
                if (!haploid_genomes[i].n)
                    rv.push(i);
            }
        return flagged_haploid_genome_queue(std::move(rv));
    }

    template <typename GenomeContainerType>
    inline std::size_t
    recycle_haploid_genome(
        GenomeContainerType &haploid_genomes,
        flagged_haploid_genome_queue &haploid_genome_recycling_bin,
        typename GenomeContainerType::value_type::mutation_container &neutral,
        typename GenomeContainerType::value_type::mutation_container &selected)
    /// \brief Return location of a new haploid_genome, recycling available memory if possible
    /// \param haploid_genomes vector of haploid_genomes
    /// \param haploid_genome_recycling_bin A flagged_haploid_genome_queue
    /// \param neutral Data for new haploid_genome's neutral variants
    /// \param selected Data for new haploid_genome's selected variants
    /// \return A location in \a haploid_genomes
    {
        // Try to recycle
        auto &ref = haploid_genome_recycling_bin.get();
        if (!ref.empty())
            {
                auto idx = ref.front();
                ref.pop();
#ifndef NDEBUG
                if (haploid_genomes[idx].n)
                    {
                        throw std::runtime_error(
                            "FWDPP DEBUG: attempting to recycle an extant "
                            "haploid_genome");
                    }
#endif
                haploid_genomes[idx].mutations.swap(neutral);
                haploid_genomes[idx].smutations.swap(selected);
                return idx;
            }
        haploid_genomes.emplace_back(0u, std::move(neutral), std::move(selected));
        return (haploid_genomes.size() - 1);
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
          MutationContainerType::value_type
         */
    template <typename MutationContainerType, class... Args>
    inline std::size_t
    recycle_mutation_helper(flagged_mutation_queue &mutation_recycling_bin,
                            MutationContainerType &mutations, Args &&... args)
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
                mutations[rv] = typename MutationContainerType::value_type(
                    std::forward<Args>(args)...);
                return rv;
            }
        mutations.emplace_back(std::forward<Args>(args)...);
        return mutations.size() - 1;
    }
} // namespace fwdpp

#endif
