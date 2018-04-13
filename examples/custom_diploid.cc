/// \include custom_diploid.cc
/// Example of a custom diploid type.

#include <sstream>
#include <stdexcept>
#include <fwdpp/type_traits.hpp>
#include <fwdpp/io/diploid.hpp>
#include <cstddef>

struct dip_with_deme
{
    using first_type = std::size_t;
    using second_type = std::size_t;
    std::size_t first, second;
    std::int32_t deme;

    dip_with_deme() : first{}, second{}, deme{} {}
    dip_with_deme(first_type f, second_type s)
        : first{ f }, second{ s }, deme{ 0 }
    {
    }
    dip_with_deme(first_type f, second_type s, std::int32_t d)
        : first{ f }, second{ s }, deme{ d }
    {
    }
};

static_assert(fwdpp::traits::is_diploid<dip_with_deme>::value,
              "Error: dip with deme is not a valid diploid type");

namespace fwdpp
{
    namespace io
    {
        template <> struct serialize_diploid<dip_with_deme>
        {
            template <typename ostreamtype>
            inline void
            operator()(ostreamtype& buffer, const dip_with_deme& dip) const
            {
                scalar_writer writer;
                writer(buffer, &dip.first);
                writer(buffer, &dip.second);
                writer(buffer, &dip.deme);
            }
        };

        template <> struct deserialize_diploid<dip_with_deme>
        {
            template <typename istreamtype>
            inline void
            operator()(istreamtype& buffer, dip_with_deme& d) const
            {
                scalar_reader reader;
                reader(buffer, &d.first);
                reader(buffer, &d.second);
                reader(buffer, &d.deme);
            }
        };
    }
}

int
main(int argc, char** argv)
{
    dip_with_deme dip(0, 34, 11);
    dip_with_deme dip2;

    std::ostringstream obuffer;
    fwdpp::io::serialize_diploid<dip_with_deme> serializer;
    serializer(obuffer, dip);

    std::istringstream ibuffer(obuffer.str());
    fwdpp::io::deserialize_diploid<dip_with_deme> deserializer;
    deserializer(ibuffer, dip2);

    if (dip.first != dip2.first || dip.second != dip2.second
        || dip.deme != dip2.deme)
        {
            throw std::runtime_error("serialization_error");
        }
}
