#ifndef FWDPP_IO_DIPLOID_HPP__
#define FWDPP_IO_DIPLOID_HPP__

#include "scalar_serialization.hpp"

namespace fwdpp
{
    namespace io
    {
        template <typename T> struct serialize_diploid
        {
            io::scalar_writer writer;
            serialize_diploid() : writer{} {}
            template <typename ostreamtype>
            inline void
            operator()(const T &dip, ostreamtype &buffer) const
            {
                writer(buffer, &dip.first);
                writer(buffer, &dip.second);
            }
        };

        template <typename T> struct deserialize_diploid
        {
            io::scalar_reader reader;
            deserialize_diploid() : reader{} {}
            template <typename istreamtype>
            inline void
            operator()(T &dip, istreamtype &buffer) const
            {
                typename T::first_type c;
                reader(buffer, &c);
                dip.first = c;
                reader(buffer, &c);
                dip.second = c;
            }
        };
    }
}

#endif
