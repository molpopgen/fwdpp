# Tree sequences

fwdpp 0.7.0 introduced support for generating "succinct tree sequences" during simulations. See \cite Kelleher2016-cb
for background on the data structures and \cite Kelleher2018-fu on their application to forward simulation.

Some readers will probably find the [msprime documentation](https://msprime.readthedocs.io/en/stable/) helpful here, as
many of the big-picture concepts are the same.  The relevant section of that documentation is titled "Tree sequence interchange".

## Example programs

I encourage you to read the [tutorials](https://tskit-dev.github.io/tutorials/) that we have written for the tskit project.  The tutorials on forward-time simulation give an overview of the logic of a simulation.  While the overall idea of "generate nodes and edges, simplify every now and then, repeat" is straightforward, there are some gotchas that tripped up the authors of \cite Kelleher2018-fu more than a few times.

1. wfts.cc implements a Wright-Fisher simulation with selection, ancestral sample recording, and ancestral sample
   metadata tracking. Command line switches also enable extensive testing of the data structures.
2. spatialts.cc implements a simplistic method of simulating discrete generations on a continuous landscape.  The
   methods used in this example are **not** what one would do in research-quality code.  Rather, the example is
   illustrative, and shows how the geographic locations of preserved samples may be tracked.

## Namespace

Tree sequence support is found in namespace fwdpp::ts. Like the rest of fwdpp, it is header-only.

## Data types

1. Nodes and edges are defined by fwdpp::ts::node and fwdpp::ts::edge, respectively.
2. Mutations are tracked on trees via fwdpp::ts::mutation_record.

Node, edge, and mutation tables are stored as vectors of the above three types.

fwdpp::ts::table_collection holds node, edge, and mutation tables.

fwdpp::ts::table_simplifier is responsible for applying the simplification algorithm of \cite Kelleher2018-fu to a
fwdpp::ts::table_collection.

## Relationship to msprime/tskit

I'll refer to *tskit* here, which means the part of msprime dealing with tree sequences.

1. Node times are measured forwards in time. See fwdpp::ts::node.
2. The mutation table is simply a vector of records tracking nodes where mutations occur on the tree and their indexes in the population.  See
   fwdpp::ts::mutation_record.  Thus, the concept of a "sites table" and a "mutation table" are merged in fwdpp.
3. The default storage format is a custom format.  Once a standalone tskit C library is available, we will likely
   provide direct access to storage in that tree file format as well.  The current binary format is required for things
   like MPI, Python pickling, etc., though, so it'll stick around.

## Iterating over trees

fwdpp::ts::tree_visitor iterates over the non-recombining segments ("marginal trees") in the pedigree.  This type
defines a call operator that handles the iteration, updating an internal variable of type fwdpp::ts::marginal_tree.

![A simplified, fully-coalesced marginal tree](images/tree.png)

Consider the tree shown above. It has 7 nodes. The parent list of the tree is represented as a linear vector
(fwdpp::ts::marginal_tree::parents), and would contain the following values:

```
4 4 5 5 6 6 fwdpp::ts::TS_NULL_NODE
```

The final value, which is the parent of the node labelled **6**, is fwdpp::ts::TS_NULL_NODE, which signifies that there
is no parent.

Iteration over parents is trivial:

```cpp
void visit_parents(const fwdpp::ts::TS_NODE_INT n, const fwdpp::ts::marginal_tree & m)
{
    auto p = n;
    while (p != fwdpp::ts::TS_NULL_NODE)
        {
            p = m.parents[p];
        }
}
```

Of course, a real-world application would have to provide bounds-checking, and would probably actually *do* something
with `p` at each iteration.

### Traversing the nodes in a tree

In general, a node table is as least as large as the number of nodes in any one tree.  We therefore need a way to access the nodes
that actually *are* in a given tree.  The class fwdpp::ts::node_iterator facilitates efficient traversal of the nodes in a fwdpp::ts::marginal_tree.  Direct interaction with fwdpp::ts:node_iterator is often unnecessary.  It will often be more convenient to use the function fwdpp::ts::process_nodes.  For example, to obtain a mapping of nodes to the number of children descending from each node:

```cpp
#include <vector>
#include <utility>
#include <fwdpp/ts/marginal_tree_functions.hpp>

std::vector<std::pair<fwdpp::ts::TS_NODE_INT, int>>
nchildren_per_node(const fwdpp::ts::marginal_tree & m)
{
	std::vector<std::pair<fwdpp::ts::TS_NODE_INT, int>> rv;
	fwdpp::ts::process_nodes(m,
				 fwdpp::ts::nodes_preorder(),
				 [&m,&rv](fwdpp::ts::TS_NODE_INT u)
				 {
					 rv.emplace_back(u, fwdpp::ts::num_children(m,u));
				 });
	return rv;
}
```

Internally, the above code constructs a fwdpp::ts::node_iterator object that conducts a preorder traversal of \a m.  The type fwdpp::ts::nodes_preorder is a dispatch tag that facilitates the dependency injection of the traversal order policy when constructing the fwdpp::ts::node_iterator.  The lambda function gets the number of children for each node via a call to fwdpp::ts::num_children (which dispatches its own work to fwdpp::ts::child_iterator).

The library provide classes and functions for iterating over the tree roots and the children of nodes:

* fwdpp::ts::root_iterator allows left-to-right traversal of all nodes in a fwdpp::ts::marginal_tree
* fwdpp::ts::num_roots returns the number of roots in a fwdpp::ts::marginal_tree
* fwdpp::ts::get_roots returns a vector of the roots of the tree.
* fwdpp::ts::process_roots allows a function to be applied to all roots of a tree.
* fwdpp::ts::child_iterator allows left-to-right or right-to-left iteration over all children of a node.
* fwdpp::ts::num_children returns the number of children descending from a node
* fwdpp::ts::get_children returns a vector of the children descending from a node
* fwdpp::ts::process_children allows a function to be applied to all children of a node.

The classes and functions listed above may be included via fwdpp/ts/marginal_tree_functions.hpp or via the individual headers included therein.

#### Node traversal order

The node traversal method is determined by a dependency injection into fwdpp::ts::node_iterator.  The dependency must publicly inherit from the abstract class fwdpp::ts::node_traversal_order.  The construction of a fwdpp::ts::node_iterator object is done via tag dispatch to a function called fwdpp::ts::node_traversal_dispatch.  For example, setting up preorder traversal works like this:

```cpp
fwdpp::ts::node_iterator mi(m, fwdpp::ts::nodes_preorder());
```

The reason to use tag dispatch here is that you may implement a traversal order that fwdpp does not provide.  To do so, you have to include fwdpp/ts/marginal_tree_functions/node_traversal_order.hpp and define your new class, a tag dispatch struct and an overload of node_traversal_dispatch **before** including fwdpp/ts/marginal_tree_functions.hpp or fwdpp/ts/marginal_tree_functions/nodes.hpp.  For a concrete example, see the implementation of fwdpp::ts::node_traversal_preorder.  Note that you do not need to put your own classes and functions into the fwdpp namespaces because argument-dependent lookup (ADL) should apply here, which is tested in the test suite.
