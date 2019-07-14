#ifndef FWDPP_TS_TABLE_TYPES_HPP
#define FWDPP_TS_TABLE_TYPES_HPP

#include <vector>
#include "node.hpp"
#include "edge.hpp"
#include "site.hpp"
#include "mutation_record.hpp"

namespace fwdpp
{
    namespace ts
    {
        /// An "edge table"
        ///  \version 0.7.0 Added to fwdpp
        using edge_vector = std::vector<edge>;
        /// A "node table"
        ///  \version 0.7.0 Added to fwdpp
        using node_vector = std::vector<node>;
        /// A "mutation table"
        ///  \version 0.7.0 Added to fwdpp
        using mutation_key_vector = std::vector<mutation_record>;
        /// \brief Site table
        /// \version 0.8.0 Added to fwdpp
        using site_vector = std::vector<site>;
    }
}

#endif
