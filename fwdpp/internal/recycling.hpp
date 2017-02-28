#ifndef FWDPP_INTERNAL_RECYCLING
#define FWDPP_INTERNAL_RECYCLING

#include <cassert>
#include <queue>
#include <cassert>
#include <type_traits>

namespace KTfwd
{
    namespace fwdpp_internal
    {

		template <class T> using recycling_bin_t = typename std::conditional<std::is_unsigned<T>::value,std::queue<T>,void>::type;

        template <typename mcount_vec>
        recycling_bin_t<typename mcount_vec::size_type>
        make_mut_queue(const mcount_vec &mcounts)
        {
            recycling_bin_t<typename mcount_vec::size_type> rv;
            const auto msize = mcounts.size();
            for (typename mcount_vec::size_type i = 0; i < msize; ++i)
                {
                    if (!mcounts[i])
                        rv.push(i);
                }
            return rv;
        }

        template <typename gvec_t>
        recycling_bin_t<typename gvec_t::size_type>
        make_gamete_queue(const gvec_t &gametes)
        {
            recycling_bin_t<typename gvec_t::size_type> rv;
            const auto gsize = gametes.size();
            for (typename gvec_t::size_type i = 0; i < gsize; ++i)
                {
                    if (!gametes[i].n)
                        rv.push(i);
                }
            return rv;
        }

        template <typename gcont_t, typename queue_t>
        inline typename queue_t::value_type
        recycle_gamete(
            gcont_t &gametes, queue_t &gamete_recycling_bin,
            typename gcont_t::value_type::mutation_container &neutral,
            typename gcont_t::value_type::mutation_container &selected)
        {
            // Try to recycle
            if (!gamete_recycling_bin.empty())
                {
                    auto idx = gamete_recycling_bin.front();
                    gamete_recycling_bin.pop();
                    assert(!gametes[idx].n);
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
        template <typename queue_t, typename mcont_t, class... Args>
        typename queue_t::value_type
        recycle_mutation_helper(queue_t &mutation_recycling_bin,
                                mcont_t &mutations, Args &&... args)
        {
            if (!mutation_recycling_bin.empty())
                {
                    auto rv = mutation_recycling_bin.front();
                    mutation_recycling_bin.pop();
                    mutations[rv] = typename mcont_t::value_type(std::forward<Args>(args)...);
                    return rv;
                }
            mutations.emplace_back(std::forward<Args>(args)...);
            return mutations.size() - 1;
        }
    }
}

#endif
