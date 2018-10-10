#ifndef FWDPP_TS_MSPRIME_ALGO_HPP
#define FWDPP_TS_MSPRIME_ALGO_HPP

/// Implements Algorithms T and L from
/// Kelleher et al., 2016, aka "the msprime
/// paper", DOI: 10.1371/journal/pcbi.1004842
/// This implementation is experimental, and
/// is intended for integration into a future
/// release of fwdpp.
/// The most likely fate is that this code
/// end up as member functions of table_collection
/// and/or table_simplifier, so that we can
/// re-use memory allocated for I, O, and pi,
/// as we will be calling these functions a LOT.
/// LICENSE: GPL3
/// Author: Kevin Thornton
/// Thanks: Jerome Kelleher

#include <vector>
#include <cstdint>
#include <type_traits>
#include "marginal_tree.hpp"
#include "iterate_marginal_trees.hpp"

namespace fwdpp
{
    namespace ts
    {
        template <typename visitor>
        void
        algorithmT(const indexed_edge_container& input_left,
                   const indexed_edge_container& output_right,
                   const std::int32_t nnodes, const double maxpos, visitor v)
        {
            marginal_tree marginal(nnodes);
            iterate_marginal_trees(marginal, input_left, output_right, maxpos,
                                   v, std::false_type(), std::false_type());
        }

        template <typename visitor>
        void
        algorithmL(const indexed_edge_container& input_left,
                   const indexed_edge_container& output_right,
                   const std::vector<std::int32_t>& sample_indexes,
                   const std::int32_t nnodes, const double maxpos, visitor v)
        {
            marginal_tree marginal(nnodes, sample_indexes);
            iterate_marginal_trees(marginal, input_left, output_right, maxpos,
                                   v, std::true_type(), std::false_type());
        }

        template <typename visitor>
        void
        algorithmL(const indexed_edge_container& input_left,
                   const indexed_edge_container& output_right,
                   const std::vector<std::int32_t>& sample_indexes,
                   const std::vector<std::int32_t>& preserved_nodes,
                   const std::int32_t nnodes, const double maxpos, visitor v)
        {
            marginal_tree marginal(nnodes, sample_indexes, preserved_nodes);
            iterate_marginal_trees(marginal, input_left, output_right, maxpos,
                                   v, std::true_type(), std::false_type());
        }

        template <typename visitor>
        void
        algorithmS(const indexed_edge_container& input_left,
                   const indexed_edge_container& output_right,
                   const std::vector<std::int32_t>& sample_indexes,
                   const std::int32_t nnodes, const double maxpos, visitor v)
        {
            marginal_tree marginal(nnodes, sample_indexes);
            iterate_marginal_trees(marginal, input_left, output_right, maxpos,
                                   v, std::false_type(), std::true_type());
        }
    } // namespace ts
} // namespace fwdpp
#endif
