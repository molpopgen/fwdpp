#include <sstream>
#include <stdexcept>
#include <fwdpp/type_traits.hpp>
#include <fwdpp/io/diploid.hpp>
#include <cstddef>

struct dip_with_fitness
{
    using first_type = std::size_t;
    using second_type = std::size_t;
    std::size_t first, second;
    double fitness;

    dip_with_fitness() : first{}, second{}, fitness{}
    {
    }
    dip_with_fitness(first_type f, second_type s) : first{f}, second{s}, fitness{0}
    {
    }
    dip_with_fitness(first_type f, second_type s, double w)
        : first{f}, second{s}, fitness{w}
    {
    }
};

static_assert(fwdpp::traits::is_diploid_v<dip_with_fitness>,
              "Error: dip with fitness is not a valid diploid type");
static_assert(fwdpp::traits::is_custom_diploid_v<dip_with_fitness>,
              "Error: dip with fitness is not a valid diploid type");

namespace fwdpp
{
    namespace io
    {
        template <> struct serialize_diploid<dip_with_fitness>
        {
            template <typename ostreamtype>
            inline void
            operator()(ostreamtype& buffer, const dip_with_fitness& dip) const
            {
                scalar_writer writer;
                writer(buffer, &dip.first);
                writer(buffer, &dip.second);
                writer(buffer, &dip.fitness);
            }
        };

        template <> struct deserialize_diploid<dip_with_fitness>
        {
            template <typename istreamtype>
            inline void
            operator()(istreamtype& buffer, dip_with_fitness& d) const
            {
                scalar_reader reader;
                reader(buffer, &d.first);
                reader(buffer, &d.second);
                reader(buffer, &d.fitness);
            }
        };
    }
}

int
main(int, char**)
{
    dip_with_fitness dip(0, 34, 11);
    dip_with_fitness dip2;

    std::ostringstream obuffer;
    fwdpp::io::serialize_diploid<dip_with_fitness> serializer;
    serializer(obuffer, dip);

    std::istringstream ibuffer(obuffer.str());
    fwdpp::io::deserialize_diploid<dip_with_fitness> deserializer;
    deserializer(ibuffer, dip2);

    if (dip.first != dip2.first || dip.second != dip2.second || dip.fitness != dip2.fitness)
        {
            throw std::runtime_error("serialization_error");
        }
}
