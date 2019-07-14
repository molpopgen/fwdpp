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
             *  \version 0.7.4 Updated value to 2
             *  \version 0.8.0 Updated value to 3
			 */
            constexpr const std::uint32_t TS_TABLES_VERSION = 3;
        }
    } // namespace ts
} // namespace fwdpp

#endif
