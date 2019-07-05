#ifndef FWDPP_TESTSUITE_SIMPLE_TABLE_COLLECTION_HPP
#define FWDPP_TESTSUITE_SIMPLE_TABLE_COLLECTION_HPP

#include <fwdpp/ts/table_collection.hpp>
#include <fwdpp/ts/tree_visitor.hpp>

class simple_table_collection
{
  private:
    fwdpp::ts::table_collection
    init_tables()
    {
        fwdpp::ts::table_collection t(1.);
        t.push_back_node(3, 0);
        t.push_back_node(3, 0);
        t.push_back_node(3, 0);
        t.push_back_node(3, 0);
        t.push_back_node(2, 0);
        t.push_back_node(1, 0);
        t.push_back_node(0, 0);
        t.push_back_edge(0, 1, 6, 5);
        t.push_back_edge(0, 1, 6, 4);
        t.push_back_edge(0, 1, 5, 2);
        t.push_back_edge(0, 1, 5, 3);
        t.push_back_edge(0, 1, 4, 1);
        t.push_back_edge(0, 1, 4, 0);
        t.sort_edges();
		t.build_indexes(); //NOTE: critical!
        return t;
    }

  public:
    fwdpp::ts::table_collection tables;
    std::vector<fwdpp::ts::TS_NODE_INT> samples;
	// NOTE: tv is auto advanced to the first tree
	// by the constructor
    fwdpp::ts::tree_visitor tv;
    explicit simple_table_collection()
        : tables(init_tables()), samples{ { 0, 1, 2, 3 } }, tv(tables, samples)
    {
        tv(std::false_type(), std::false_type());
    }
};

#endif
