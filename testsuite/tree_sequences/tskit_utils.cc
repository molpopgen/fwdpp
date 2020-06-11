#include "tskit_utils.hpp"
#include <stdexcept>
#include <sstream>

void
handle_tskit_return_code(int code)
{
    if (code < 0)
        {
            std::ostringstream o;
            o << tsk_strerror(code);
            throw std::runtime_error(o.str());
        }
}

tsk_treeseq_wrapper::tsk_treeseq_wrapper(tsk_table_collection_t* tables) : treeseq{}
{
    auto rv = tsk_treeseq_init(&treeseq, tables, 0);
    handle_tskit_return_code(rv);
}

tsk_treeseq_wrapper::~tsk_treeseq_wrapper()
{
    auto rv = tsk_treeseq_free(&treeseq);
    handle_tskit_return_code(rv);
}

tsk_treeseq_t*
tsk_treeseq_wrapper::get()
{
    return &treeseq;
};

const tsk_treeseq_t*
tsk_treeseq_wrapper::get() const
{
    return &treeseq;
};

table_collection_ptr
make_table_collection_ptr(double sequence_length)
{
    table_collection_ptr rv(new tsk_table_collection_t(),
                            [](tsk_table_collection_t* tables) {
                                tsk_table_collection_free(tables);
                                delete tables;
                            });
    int err = tsk_table_collection_init(rv.get(), 0);
    handle_tskit_return_code(err);
    rv->sequence_length = sequence_length;
    if (err != 0)
        {
            throw std::runtime_error("could not initialize tsk_table_collection");
        }
    return rv;
}

table_collection_ptr
dump_table_collection_to_tskit(const fwdpp::ts::std_table_collection& tables,
                               double forward_time, std::vector<int>& samples)
{
    if (samples.size() != tables.num_nodes())
        {
            throw std::runtime_error("samples.size() != tables.num_nodes()");
        }
    auto tskit_tables = make_table_collection_ptr(tables.genome_length());
    unsigned i = 0;
    for (auto& n : tables.nodes)
        {
            // Convert time from forwards to backwards
            // and label the last generation (first 2N nodes)
            // as samples.
            int rv = tsk_node_table_add_row(&tskit_tables->nodes,
                                            (samples[i] != 0) ? TSK_NODE_IS_SAMPLE
                                                              : 0,         // flag
                                            -1. * (n.time - forward_time), // time
                                            TSK_NULL,                      // population
                                            TSK_NULL,                      // individual
                                            nullptr,                       // metadata
                                            0); // metadata length
            handle_tskit_return_code(rv);
            ++i;
        }
    for (auto& e : tables.edges)
        {
            auto rv = tsk_edge_table_add_row(&tskit_tables->edges, e.left, e.right,
                                             e.parent, e.child, nullptr, 0);
            handle_tskit_return_code(rv);
        }
    auto rv = tsk_table_collection_build_index(tskit_tables.get(), 0);
    handle_tskit_return_code(rv);
    return tskit_tables;
}

