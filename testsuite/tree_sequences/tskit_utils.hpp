#ifndef FWDPP_TESTSUITE_TREE_SEQUENCES_TSKIT_UTILS_HPP
#define FWDPP_TESTSUITE_TREE_SEQUENCES_TSKIT_UTILS_HPP

#include <cstdint>
#include <cstdio>
#include <tskit.h>
#include <memory>
#include <functional>
#include <fwdpp/ts/std_table_collection.hpp>

using table_collection_ptr
    = std::unique_ptr<tsk_table_collection_t,
                      std::function<void(tsk_table_collection_t*)>>;

table_collection_ptr make_table_collection_ptr(double sequence_length);

table_collection_ptr
dump_table_collection_to_tskit(const fwdpp::ts::std_table_collection& tables,
                               double forward_time, std::vector<int> & samples);

struct tsk_treeseq_wrapper
{
	tsk_treeseq_t treeseq;
	explicit tsk_treeseq_wrapper(tsk_table_collection_t * tables);
	~tsk_treeseq_wrapper();

	tsk_treeseq_t * get();
	const tsk_treeseq_t * get() const;
};

void handle_tskit_return_code(int code);
#endif
