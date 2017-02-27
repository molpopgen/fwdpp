#ifndef FWDPP_UNIT_TESTING_CUSTOM_DIP_HPP
#define FWDPP_UNIT_TESTING_CUSTOM_DIP_HPP

#include <iosfwd>
#include <limits>
#include <fwdpp/type_traits.hpp>

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

static_assert(KTfwd::traits::is_custom_diploid_t<custom_diploid_testing_t>::value,
		"custom_diploid_testing_t must pass as a custom diploid.");

struct diploid_writer
{
    using result_type = void;
    template <typename itr, typename streamtype>
    inline result_type
    operator()(itr i, streamtype &o) const
    {
        o.write(reinterpret_cast<const char *>(&i.i), sizeof(unsigned));
    }
};

struct diploid_reader
{
    using result_type = void;
    template <typename itr, typename streamtype>
    inline result_type
    operator()(itr i, streamtype &in) const
    {
        in.read(reinterpret_cast<char *>(&i.i), sizeof(unsigned));
    }
};

#endif
