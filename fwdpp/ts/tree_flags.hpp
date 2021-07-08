#pragma once

#include <cstdint>

namespace fwdpp
{
    namespace ts
    {
        namespace tree_flags
        {
            static constexpr std::uint32_t NONE = 0;
            static constexpr std::uint32_t TRACK_SAMPLES = 1 << 0;
        }
    }
}
