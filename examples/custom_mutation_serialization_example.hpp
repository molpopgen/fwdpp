#pragma once

#include "custom_mutation_example.hpp"
#include <fwdpp/io/mutation.hpp>

namespace fwdpp
{
    namespace io
    {
        template <> struct serialize_mutation<mutation>
        {
            template <typename streamtype>
            inline void
            operator()(streamtype &buffer, const mutation &m) const
            {
                io::scalar_writer writer;
                writer(buffer, &m.pos);
                writer(buffer, &m.s);
                writer(buffer, &m.h);
            }
        };

        template <> struct deserialize_mutation<mutation>
        {
            template <typename streamtype>
            inline mutation
            operator()(streamtype &buffer) const
            {
                double pos, s, h;
                io::scalar_reader reader;
                reader(buffer, &pos);
                reader(buffer, &s);
                reader(buffer, &h);

                return mutation(pos, s, h);
            }
        };

    }
}
