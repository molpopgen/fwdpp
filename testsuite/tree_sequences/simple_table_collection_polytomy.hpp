#ifndef FWDPP_TESTSUITE_SIMPLE_TABLE_COLLECTION_POLYTOMY_HPP
#define FWDPP_TESTSUITE_SIMPLE_TABLE_COLLECTION_POLYTOMY_HPP

#include <fwdpp/ts/table_collection.hpp>
#include <fwdpp/ts/tree_visitor.hpp>

class simple_table_collection_polytomy
//         7
//      --------
//  	|      |
//  	|      6
//  	5     ---
//	   -----  | |
//	   | | |  | |
//	   0 1 2  3 4
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
        t.push_back_node(3, 0); // Node 4

        t.push_back_node(2, 0); // Node 5
        t.push_back_node(1, 0); // Node 6
        t.push_back_node(0, 0); // Node 7
        t.push_back_edge(0, 1, 7, 5);
        t.push_back_edge(0, 1, 7, 6);
        t.push_back_edge(0, 1, 6, 4);
        t.push_back_edge(0, 1, 6, 3);
        t.push_back_edge(0, 1, 5, 0);
        t.push_back_edge(0, 1, 5, 1);
        t.push_back_edge(0, 1, 5, 2);
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
    double total_time;
    explicit simple_table_collection_polytomy()
        : tables(init_tables()), samples{ { 0, 1, 2, 3, 4 } },
          tv(tables, samples, fwdpp::ts::update_samples_list(true)),
          total_time(10)
    {
        tv();
    }

    void
    reset_visitor(bool update_samples_list)
    {
        tv = fwdpp::ts::tree_visitor(
            tables, samples,
            fwdpp::ts::update_samples_list(update_samples_list));
        tv();
    }

    void
    reset_visitor_and_samples(
        const std::vector<fwdpp::ts::TS_NODE_INT>& new_samples_list,
        bool update_samples_list)
    {
        samples = new_samples_list;
        tv = fwdpp::ts::tree_visitor(
            tables, new_samples_list,
            fwdpp::ts::update_samples_list(update_samples_list));
        tv();
    }
};

#endif

