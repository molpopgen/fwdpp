#ifndef FWDPP_TS_SIMPLIFY_TABLES_OUTPUT_HPP
#define FWDPP_TS_SIMPLIFY_TABLES_OUTPUT_HPP

#include <vector>
#include <type_traits>
#include <fwdpp/ts/definitions.hpp>

namespace fwdpp
{
    namespace ts
    {
        template <typename NodeIDMap, typename MutationIndexVector>
        struct simplify_tables_output_t
        {
            static_assert(std::is_integral<typename NodeIDMap::value_type>::value,
                          "NodeIDMap::value_type must be integral");
            static_assert(std::is_signed<typename NodeIDMap::value_type>::value,
                          "NodeIDMap::value_type must be signed");
            static_assert(sizeof(typename NodeIDMap::value_type)
                              >= sizeof(table_index_t),
                          "sizeof(NodeIDMap::value_type) must be >= "
                          "sizeof(fwdpp::ts::table_index_t)");
            NodeIDMap idmap;
            MutationIndexVector preserved_mutations;

            simplify_tables_output_t() : idmap{}, preserved_mutations{}
            {
            }

            void
            clear()
            {
                idmap.clear();
                preserved_mutations.clear();
            }
        };

        using simplify_tables_output
            = simplify_tables_output_t<std::vector<table_index_t>,
                                       std::vector<std::size_t>>;
    }
}

#endif
