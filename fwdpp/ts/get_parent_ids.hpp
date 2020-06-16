#ifndef FWDPP_TS_GET_PARENT_IDS_HPP
#define FWDPP_TS_GET_PARENT_IDS_HPP

#include <utility>
#include "definitions.hpp"

namespace fwdpp
{
    namespace ts
    {
        inline std::pair<table_index_t, table_index_t>
        get_parent_ids(const table_index_t first_parental_index,
                       const table_index_t parent, const int did_swap)
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
            return std::make_pair(
                first_parental_index + 2 * static_cast<table_index_t>(parent)
                    + did_swap,
                first_parental_index + 2 * static_cast<table_index_t>(parent)
                    + !did_swap);
        }
    } // namespace ts
} // namespace fwdpp

#endif
