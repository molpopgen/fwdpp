#include <stack>
#include <vector>
#include <memory>
#include <fwdpp/ts/marginal_tree_functions/node_traversal_order.hpp>
#include <fwdpp/ts/marginal_tree_functions/children.hpp>

class node_traversal_preorder_test : public fwdpp::ts::node_traversal_order
/// This is an exact copy of fwdpp::ts::node_traversal_preorder,
/// moved into gloal namespace to mimic client-side code using ADL +
/// tag dispatch
{
  private:
    using node_stack = std::stack<fwdpp::ts::TS_NODE_INT,
                                  std::vector<fwdpp::ts::TS_NODE_INT>>;
    node_stack nstack;
    fwdpp::ts::TS_NODE_INT current_node;

  public:
    explicit node_traversal_preorder_test(fwdpp::ts::TS_NODE_INT u)
        : node_traversal_order(), nstack{ { u } }, current_node{
              fwdpp::ts::TS_NULL_NODE
          }
    {
    }

    fwdpp::ts::TS_NODE_INT
    operator()(const fwdpp::ts::marginal_tree& m) final
    {
        if (!nstack.empty())
            {
                current_node = nstack.top();
                nstack.pop();
                if (num_children(m, current_node) != 0)
                    {
                        process_children(m, current_node, false,
                                         [this](fwdpp::ts::TS_NODE_INT x) {
                                             nstack.push(x);
                                         });
                    }
                return current_node;
            }
        return fwdpp::ts::TS_NULL_NODE;
    }

    void
    initialize(fwdpp::ts::TS_NODE_INT root) final
    {
        nstack.push(root);
    }
};

struct nodes_preorder_test
{
};

inline std::unique_ptr<fwdpp::ts::node_traversal_order>
node_traversal_dispatch(fwdpp::ts::TS_NODE_INT root, nodes_preorder_test)
{
    return std::unique_ptr<fwdpp::ts::node_traversal_order>(
        new node_traversal_preorder_test(root));
}
