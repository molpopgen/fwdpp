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
			 */
            constexpr const std::uint32_t TS_TABLES_VERSION = 2;
        }
    } // namespace ts
} // namespace fwdpp

#endif
