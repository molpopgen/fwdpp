#include <stack>
#include <vector>
#include <memory>
#include <fwdpp/ts/marginal_tree_functions/node_traversal_order.hpp>
#include <fwdpp/ts/marginal_tree_functions/children.hpp>

class node_traversal_preorder_test : public fwdpp::ts::node_traversal_order<std::int32_t>
/// This is an exact copy of fwdpp::ts::node_traversal_preorder,
/// moved into gloal namespace to mimic client-side code using ADL +
/// tag dispatch
{
  private:
    using node_stack = std::stack<std::int32_t, std::vector<std::int32_t>>;
    node_stack nstack;
    std::int32_t current_node;

  public:
    explicit node_traversal_preorder_test(std::int32_t u)
        : node_traversal_order(), nstack{{u}}, current_node{-1}
    {
    }

    std::int32_t
    operator()(const fwdpp::ts::marginal_tree<std::int32_t>& m) final
    {
        if (!nstack.empty())
            {
                current_node = nstack.top();
                nstack.pop();
                if (num_children(m, current_node) != 0)
                    {
                        process_children(m, current_node, false,
                                         [this](std::int32_t x) { nstack.push(x); });
                    }
                return current_node;
            }
        return -1;
    }

    void
    initialize(std::int32_t root) final
    {
        nstack.push(root);
    }
};

struct nodes_preorder_test
{
};

inline std::unique_ptr<fwdpp::ts::node_traversal_order<std::int32_t>>
node_traversal_dispatch(std::int32_t root, nodes_preorder_test)
{
    return std::unique_ptr<fwdpp::ts::node_traversal_order<std::int32_t>>(
        new node_traversal_preorder_test(root));
}
