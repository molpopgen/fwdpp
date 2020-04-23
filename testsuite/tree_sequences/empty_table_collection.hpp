#ifndef FWDPP_TESTSUITE_EMPTY_TABLE_COLLECTION_HPP
#define FWDPP_TESTSUITE_EMPTY_TABLE_COLLECTION_HPP

#include <fwdpp/ts/std_table_collection.hpp>

struct empty_table_collection
{
    fwdpp::ts::std_table_collection tables;
    std::vector<fwdpp::ts::TS_NODE_INT> empty_samples;

    empty_table_collection() : tables(1.), empty_samples{} {}
};

#endif
