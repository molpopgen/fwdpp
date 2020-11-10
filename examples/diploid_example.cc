#include <utility>
#include <cstdint>
#include <fwdpp/type_traits.hpp>

using diploid_size_t = std::pair<std::size_t, std::size_t>;
using diploid_int32 = std::pair<std::int32_t, std::int32_t>;

int
main(int, char **)
{
    static_assert(fwdpp::traits::is_diploid_v<diploid_size_t>,
                  "diploid_size_t is not a diploid!");
    static_assert(!fwdpp::traits::is_custom_diploid_v<diploid_size_t>,
                  "diploid_size_t is a custom diploid!");
    static_assert(fwdpp::traits::is_diploid_v<diploid_int32>,
                  "diploid_int32 is not a diploid!");
    static_assert(!fwdpp::traits::is_custom_diploid_v<diploid_int32>,
                  "diploid_size_t is a custom diploid!");
}
