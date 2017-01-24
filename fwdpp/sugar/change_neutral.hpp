/*!
  \file fwdpp/sugar/change_neutral.hpp

  \brief Facilitate changing "selection coefficients" of mutations
 */
#ifndef FWDPP_SUGAR_CHANGE_NEUTRAL_HPP
#define FWDPP_SUGAR_CHANGE_NEUTRAL_HPP

#include <algorithm>
#include <cassert>
#include <exception>
#include <fwdpp/debug.hpp>

namespace KTfwd
{
    namespace sugar
    {
        template <typename mcont_t, typename mut_key_cont_t>
        void
        change_neutral_details(const mcont_t &mutations, const double pos,
                               const std::size_t mindex, mut_key_cont_t &a,
                               mut_key_cont_t &b)
        /*!
          Implementation details of change_neutral.

          \note Do not call directly.
         */
        {
            /*
              Ask if gamete has this mutation.  We do a linear
              search here in case of finite-site schemes and we
              also do not expect this function to be called a whole lot.
            */
            auto i = std::find(std::begin(a), std::end(a), mindex);
            // If it exists and is mindex, erase it
            if (i == std::end(a))
                return;

            a.erase(i);

            // Add mutation key into b
            b.insert(
                std::upper_bound(
                    std::begin(b), std::end(b), pos,
                    [&mutations](const double &p,
                                 const typename mut_key_cont_t::value_type &vt)

                    { return p < mutations[vt].pos; }),
                mindex);
        }
    }

    template <typename poptype>
    void
    change_neutral(poptype &p, const std::size_t mindex)
    /*!
      \brief Change the value of mutation_base::neutral

      \param p A population
      \param mindex The key of the mutation you wish to change.

      This function changes mutation_base::neutral to !mutation_base::neutral
      at position
      mindex and updates the storage of this mutation in all gametes.

      \note This function does not change any other member data at
      p.mutations[mindex]!!  Thus, if you change a mutation
      from "neutral" to "selected", you must manually change the relevant
      member data to reflect its
      new effect size.  For the simplest use case of making a mutation no
      longer subject to selection,
      no additional changes are needed.  Rather, it is sufficient that
      mutation_base::neutral == true.

      \throw std::out_of_range if \a mindex is out of range.

      \ingroup sugar
    */
    {
        if (mindex >= p.mutations.size())
            throw std::out_of_range("mindex >= p.mutations.size()");

        bool is_neutral = p.mutations[mindex].neutral;

        // Change neutral flag
        p.mutations[mindex].neutral = !p.mutations[mindex].neutral;
        const auto pos = p.mutations[mindex].pos;

        for (auto &g : p.gametes)
            {
                if (g.n)
                    {
                        if (is_neutral)
                            {
                                sugar::change_neutral_details(
                                    p.mutations, pos, mindex, g.mutations,
                                    g.smutations);
                            }
                        else
                            {
                                sugar::change_neutral_details(
                                    p.mutations, pos, mindex, g.smutations,
                                    g.mutations);
                            }
                        assert(gamete_data_sane(g, p.mutations, p.mcounts));
                    }
            }
    }
}

#endif
