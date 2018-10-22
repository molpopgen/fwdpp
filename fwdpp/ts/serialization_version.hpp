#ifndef FWDPP_TS_SERIALIZATION_VERSION_HPP
#define FWDPP_TS_SERIALIZATION_VERSION_HPP

#include <cstdint>
namespace fwdpp
{
    namespace ts
    {
        namespace io
        {
			/*! \brief Current version number of binary formats
			 *  \version 0.7.0 Added to library
			 */
            constexpr const std::uint32_t TS_TABLES_VERSION = 1;
        }
    } // namespace ts
} // namespace fwdpp

#endif
