#ifndef FWDPP_IO_MUTATION_HPP__
#define FWDPP_IO_MUTATION_HPP__

#include <stdexcept>

namespace fwdpp
{
    namespace io
    {
        template <typename T> struct serialize_mutation
        {
            template <typename istreamtype>
            inline void
            operator()(const T &, istreamtype &) const
            {
                throw std::runtime_error(
                    "serializtion not implemented for this mutation type");
            }
        };

        template <typename T> struct deserialize_mutation
        {
            template <typename ostreamtype>
            inline T
            operator()(ostreamtype &) const
            {
                throw std::runtime_error(
                    "deserializtion not implemented for this mutation type");
            }
        };
    }
}

#endif
