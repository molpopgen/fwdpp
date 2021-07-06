#ifndef FWDPP_TS_GET_PARENT_IDS_HPP
#define FWDPP_TS_GET_PARENT_IDS_HPP

#include <utility>
#include <type_traits>

namespace fwdpp
{
    namespace ts
    {
        template <typename SignedInteger>
        inline std::pair<SignedInteger, SignedInteger>
        get_parent_ids(const SignedInteger first_parental_index,
                       const SignedInteger parent, const int did_swap)
        /*! 
         * Convert the index of a parent into the two node IDs.
         *
         * \param first_parental_index First index of possible parents for offspring
         * \param parent Index of parent in the population
         * \param did_swap This is the "Mendel" step.
         *
         * \note This function is trivial in implementation and is easily replaced
         * if a more complex set of rules are needed
         */
        {
            static_assert(std::is_integral<SignedInteger>::value,
                          "SignedInteger must be an integral type");
            static_assert(std::is_signed<SignedInteger>::value,
                          "SignedInteger must be a signed type");
            return std::make_pair(
                first_parental_index + 2 * static_cast<SignedInteger>(parent) + did_swap,
                first_parental_index + 2 * static_cast<SignedInteger>(parent)
                    + !did_swap);
        }
    } // namespace ts
} // namespace fwdpp

#endif
