#ifndef FWDPP_SUGAR_DEMOGRAPHY_HPP
#define FWDPP_SUGAR_DEMOGRAPHY_HPP

#include <type_traits>
#include <fwdpp/demography.hpp>
#include <fwdpp/sugar/metapop.hpp>

namespace KTfwd
{
    /*!
      Ensure that KTfwd::sugar::metapop::Ns contain sizes of all
      elements in KTfwd::sugar::metapop::diploids

      \param mpop An object of type KTfwd::metapop or KTfwd::sugar::metapop
    */
    template <typename mpoptype>
    void
    update_Ns(mpoptype &mpop)
    {
        static_assert(std::is_same<typename mpoptype::popmodel_t,
                                   sugar::METAPOP_TAG>::value,
                      "mpoptype must be an object of type KTfwd::metapop or "
                      "KTfwd::sugar::metapop");
        mpop.Ns.clear();
        for (const auto &dip : mpop.diploids)
            {
                mpop.Ns.push_back(dip.size());
            }
    }

    /*! \brief Copy a deme
      Make an exact copy of the i-th population

      \param mpop An object of type KTfwd::metapop or KTfwd::sugar::metapop
      \param i The deme to copy.

      \note Implemented in terms of KTfwd::copy_deme. See documentation of that
      function for details.
    */
    template <typename mpoptype>
    int
    copy_pop(mpoptype &mpop, const size_t i)
    {
        static_assert(std::is_same<typename mpoptype::popmodel_t,
                                   sugar::METAPOP_TAG>::value,
                      "mpoptype must be an object of type KTfwd::metapop or "
                      "KTfwd::sugar::metapop");
        auto rv = copy_deme(mpop.mutations, mpop.mcounts, mpop.gametes,
                            mpop.diploids, i);
        if (rv)
            return rv;   // there was an error
        update_Ns(mpop); // update deme sizes
        return rv;
    }

    /*! \brief Merge two demes

      Merge two demes into one

      \param mpop An object of type KTfwd::metapop or KTfwd::sugar::metapop
      \param i One deme to merge
      \param j The other deme to merge

      \note Implemented in terms of KTfwd::merge_demes. See documentation of
      that function for details.
    */
    template <typename mpoptype>
    int
    merge_pops(mpoptype &mpop, const size_t i, const size_t j)
    {
        static_assert(std::is_same<typename mpoptype::popmodel_t,
                                   sugar::METAPOP_TAG>::value,
                      "mpoptype must be an object of type KTfwd::metapop or "
                      "KTfwd::sugar::metapop");
        auto rv = merge_demes(mpop.diploids, i, j);
        if (rv)
            return rv;   // there was an error
        update_Ns(mpop); // update deme sizes
        return rv;
    }

    /*! \brief Delete a deme

      Delete a deme and remove it from the metapopulation

      \param mpop An object of type KTfwd::metapop or KTfwd::sugar::metapop
      \param i The deme to remove

      \note Implemented in terms of KTfwd::remove_deme. See documentation of
      that function for details.
    */
    template <typename mpoptype>
    int
    remove_pop(mpoptype &mpop, const size_t i)
    {
        static_assert(std::is_same<typename mpoptype::popmodel_t,
                                   sugar::METAPOP_TAG>::value,
                      "mpoptype must be an object of type KTfwd::metapop or "
                      "KTfwd::sugar::metapop");
        auto rv = remove_deme(mpop.mutations, mpop.mcounts, mpop.gametes,
                              mpop.diploids, i);
        if (rv)
            return rv;   // there was an error
        update_Ns(mpop); // update deme sizes
        return rv;
    }

    /*! \brief Swap two demes

      Swap two sub-populations

      \param mpop An object of type KTfwd::metapop or KTfwd::sugar::metapop
      \param i One deme to swap
      \param j The other deme to swap

      \note Implemented in terms of KTfwd::swap_demes. See documentation of
      that function for details.
    */
    template <typename mpoptype>
    int
    swap_pops(mpoptype &mpop, const size_t i, const size_t j)
    {
        static_assert(std::is_same<typename mpoptype::popmodel_t,
                                   sugar::METAPOP_TAG>::value,
                      "mpoptype must be an object of type KTfwd::metapop or "
                      "KTfwd::sugar::metapop");
        auto rv = swap_demes(mpop.diploids, i, j);
        if (rv)
            return rv;                     // error
        std::swap(mpop.Ns[i], mpop.Ns[j]); // Faster than a call to update_Ns
        return rv;
    }

    /*! \brief "Bud" off a new sub-population

      \param r Pointer to a const gsl_rng object
      \param mpop An object of type KTfwd::metapop or KTfwd::sugar::metapop
      \param i The parental deme
      \param N_new The size of the new deme
      \param replacement Sample parents from \c i with or without replacement

      \note Implemented in terms of KTfwd::split_deme. See documentation of
      that function for details.
     */
    template <typename mpoptype>
    int
    split_pop(const gsl_rng *r, mpoptype &mpop, const size_t i,
              const uint_t N_new, const bool replacement = false)
    {
        static_assert(std::is_same<typename mpoptype::popmodel_t,
                                   sugar::METAPOP_TAG>::value,
                      "mpoptype must be an object of type KTfwd::metapop or "
                      "KTfwd::sugar::metapop");
        assert(i < mpop.diploids.size());
        auto rv = split_deme(r, mpop.mutations, mpop.mcounts, mpop.gametes,
                             mpop.diploids, i, N_new, replacement);
        if (rv)
            return rv;
        update_Ns(mpop);
        return rv;
    }

    /*! \brief Create an admixed population

    \param r Pointer to a const gsl_rng object
    \param mpop An object of type KTfwd::metapop or KTfwd::sugar::metapop
    \param i One parental deme
    \param j The other parental deme
    \param pi The fraction of parents in the new deme originating from deme \c
    i
    \param N_new The size of the new deme
    \param replacement Sample parents from \c i with or without replacement

    \note Implemented in terms of KTfwd::admix_demes. See documentation of that
    function for details.
   */
    template <typename mpoptype>
    int
    admix_pops(const gsl_rng *r, mpoptype &mpop, const size_t i,
               const size_t j, const double pi, const uint_t N_new,
               const bool replacement = false)
    {
        static_assert(std::is_same<typename mpoptype::popmodel_t,
                                   sugar::METAPOP_TAG>::value,
                      "mpoptype must be an object of type KTfwd::metapop or "
                      "KTfwd::sugar::metapop");
        auto rv = admix_demes(r, mpop.mutations, mpop.mcounts, mpop.gametes,
                              mpop.diploids, i, j, pi, N_new, replacement);
        if (rv)
            return rv;
        update_Ns(mpop);
        return rv;
    }
}

#endif
