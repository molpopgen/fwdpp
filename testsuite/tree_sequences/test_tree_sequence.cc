#include <cstdint>
#include <iostream>
#include <boost/test/unit_test.hpp>
#include <fwdpp/ts/types/tree_sequence.hpp>

using tree_sequence = fwdpp::ts::types::tree_sequence<std::int32_t>;
using tree_sequence64 = fwdpp::ts::types::tree_sequence<std::int64_t>;

//explicit instantiations
template class fwdpp::ts::types::tree_sequence<std::int32_t>;
template class fwdpp::ts::types::tree_sequence<std::int64_t>;
template class fwdpp::ts::types::tree_iterator<std::int32_t>;
template class fwdpp::ts::types::tree_iterator<std::int64_t>;

