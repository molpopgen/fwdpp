#ifndef FWDPP_TS_SIMPLIFY_TABLES_OUTPUT_HPP
#define FWDPP_TS_SIMPLIFY_TABLES_OUTPUT_HPP

#include <vector>
#include <type_traits>
#include <fwdpp/ts/definitions.hpp>

namespace fwdpp
{
    namespace ts
    {
        template <typename SignedInteger, typename UnsignedInteger>
        struct simplify_tables_output_t
        {
            static_assert(std::is_integral<SignedInteger>::value,
                          "SignedInteger must be integral");
            static_assert(std::is_signed<SignedInteger>::value,
                          "SignedInteger must be signed");
            static_assert(std::is_integral<UnsignedInteger>::value,
                          "UnsignedInteger must be integral");
            static_assert(!std::is_signed<UnsignedInteger>::value,
                          "UnsignedInteger must be signed");
            std::vector<SignedInteger> idmap;
            std::vector<UnsignedInteger> preserved_mutations;

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

        template <typename SignedInteger>
        using simplify_tables_output
            = simplify_tables_output_t<SignedInteger, std::size_t>;
    }
}

#endif
