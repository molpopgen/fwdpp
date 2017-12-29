#ifndef FWDPP_UNIT_TESTING_CUSTOM_DIP_HPP
#define FWDPP_UNIT_TESTING_CUSTOM_DIP_HPP

#include <iosfwd>
#include <limits>
#include <fwdpp/type_traits.hpp>
#include <fwdpp/IO.hpp>
// Custom diploid type.
struct custom_diploid_testing_t
/*!
  Semantically identical to standard diploid type, which is
  pair<size_t,size_t>,
  but this forces fwdpp to assume that it is a custom type, which we exploit
  for
  testing.
 */
{
    using first_type = std::size_t;
    using second_type = std::size_t;
    first_type first;
    second_type second;
    unsigned i;
    custom_diploid_testing_t()
        : first(first_type()), second(second_type()),
          i(std::numeric_limits<unsigned>::max())
    {
    }
    custom_diploid_testing_t(first_type g1, first_type g2)
        : first(g1), second(g2), i(std::numeric_limits<unsigned>::max())
    {
    }
    bool
    operator==(const custom_diploid_testing_t &rhs) const
    {
        return this->first == rhs.first && this->second == rhs.second
               && this->i == rhs.i;
    }
};

static_assert(
    fwdpp::traits::is_custom_diploid_t<custom_diploid_testing_t>::value,
    "custom_diploid_testing_t must pass as a custom diploid.");

// specialization for ADL
namespace fwdpp
{
    template <>
    inline void
    serialize_diploid<std::stringstream, custom_diploid_testing_t>(
        const custom_diploid_testing_t &dip, std::stringstream &buffer)
    {
        fwdpp::fwdpp_internal::scalar_writer w;
        w(buffer, &dip.first);
        w(buffer, &dip.second);
        w(buffer, &dip.i);
    }

    template <>
    inline void
    deserialize_diploid<std::stringstream, custom_diploid_testing_t>(
        custom_diploid_testing_t &dip, std::stringstream &buffer)
    {
        fwdpp::fwdpp_internal::scalar_reader r;
        r(buffer, &dip.first);
        r(buffer, &dip.second);
        r(buffer, &dip.i);
    }
}

#endif
