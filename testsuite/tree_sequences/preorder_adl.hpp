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
    using node_stack = std::stack<fwdpp::ts::table_index_t,
                                  std::vector<fwdpp::ts::table_index_t>>;
    node_stack nstack;
    fwdpp::ts::table_index_t current_node;

  public:
    explicit node_traversal_preorder_test(fwdpp::ts::table_index_t u)
        : node_traversal_order(), nstack{ { u } }, current_node{
              fwdpp::ts::NULL_INDEX
          }
    {
    }

    fwdpp::ts::table_index_t
    operator()(const fwdpp::ts::marginal_tree& m) final
    {
        if (!nstack.empty())
            {
                current_node = nstack.top();
                nstack.pop();
                if (num_children(m, current_node) != 0)
                    {
                        process_children(m, current_node, false,
                                         [this](fwdpp::ts::table_index_t x) {
                                             nstack.push(x);
                                         });
                    }
                return current_node;
            }
        return fwdpp::ts::NULL_INDEX;
    }

    void
    initialize(fwdpp::ts::table_index_t root) final
    {
        nstack.push(root);
    }
};

struct nodes_preorder_test
{
};

inline std::unique_ptr<fwdpp::ts::node_traversal_order>
node_traversal_dispatch(fwdpp::ts::table_index_t root, nodes_preorder_test)
{
    return std::unique_ptr<fwdpp::ts::node_traversal_order>(
        new node_traversal_preorder_test(root));
}
