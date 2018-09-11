#ifndef FWDPP_INTERNAL_DEBUG_DETAILS_HPP
#define FWDPP_INTERNAL_DEBUG_DETAILS_HPP

#include <algorithm>
#include <stdexcept>

namespace fwdpp
{
    namespace debug
    {
        namespace detail
        {
#ifndef NDEBUG
            template <typename mcont_t, typename iterator>
            void
            validate_mutation_key_ranges(const mcont_t &mutations,
                                         const iterator beg,
                                         const iterator end)
            {
                if (std::any_of(beg, end, [&mutations](const std::size_t m) {
                        return m >= mutations.size();
                    }))
                    {
                        throw std::runtime_error("FWDPP DEBUG: mutation key "
                                                 "out of range");
                    }
            }
#else
            template <typename mcont_t, typename iterator>
            void
            validate_mutation_key_ranges(const mcont_t& /*mutations*/,
                                         const iterator /*beg*/,
                                         const iterator /*end*/)
            {
            }
#endif
        }
    }
}

#endif
