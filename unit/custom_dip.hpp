#ifndef FWDPP_UNIT_TESTING_CUSTOM_DIP_HPP
#define FWDPP_UNIT_TESTING_CUSTOM_DIP_HPP

#include <fwdpp/tags/diploid_tags.hpp>
#include <iostream>

//Custom diploid type.
struct diploid_t : public KTfwd::tags::custom_diploid_t
{
  using first_type = std::size_t;
  using second_type = std::size_t;
  first_type first;
  second_type second;
  unsigned i;
  diploid_t() : first(first_type()),second(second_type()),i(std::numeric_limits<unsigned>::max()) {}
  diploid_t(first_type g1, first_type g2) : first(g1),second(g2),i(std::numeric_limits<unsigned>::max()) {}
  bool operator==(const diploid_t & rhs) const
  {
    return this->first == rhs.first &&
    this->second == rhs.second &&
    this->i == rhs.i;
  }
};

struct diploid_writer
{
  using result_type = void;
  template<typename itr,typename streamtype>
  inline result_type operator()( itr i, streamtype & o ) const
  {
    o.write( reinterpret_cast<const char *>(&i.i),sizeof(unsigned) );
  }
};

struct diploid_reader
{
  using result_type = void;
  template<typename itr,typename streamtype>
  inline result_type operator()( itr i, streamtype & in ) const
  {
    in.read( reinterpret_cast<char *>(&i.i),sizeof(unsigned) );
  }
};

#endif
