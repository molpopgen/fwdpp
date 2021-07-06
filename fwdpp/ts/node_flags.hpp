#pragma once

#include <cstdint>

namespace fwdpp
{
    namespace ts
    {
        namespace node_flags
        {
            static constexpr std::uint32_t NONE = 0;
            static constexpr std::uint32_t IS_SAMPLE = 1 << 0;
            static constexpr std::uint32_t IS_ALIVE = 1 << 1;
        }
    }
}

