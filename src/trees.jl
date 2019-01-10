export data

"""
    data(node)

Retrieve the data aka payload associated with the given node.
"""
function data end


function root end
function children end

# import AbstractTrees: print_tree, children
# export print_tree, children
#
# export insertchild, data
#
# function insertchild end
# function data end
#
#
# function depthfirst(f, t, level = 1)
#     f(t, level)
#     for c in children(t)
#         depthfirst(f, c, level+1)
#     end
# end

function depthfirst(f, tree, node=root(tree), level=1)
    f(tree, node, level)
    for c in FMMTrees.children(tree, node)
        depthfirst(f, tree, c, level+1)
    end
end

function print_tree(tree, node=root(tree))
    depthfirst(tree) do tree, node, level
        print("-"^(level-1))
        print(data(tree, node))
        println()
    end
end

"""
    update!(tree, node, data, router!, updater!)

Algorithm to update or add nodes of the tree. `router!` and `updater!` are
user supplied functions:

    router!(tree, node)

Returns the next candidate `node` until the node for insertion is reaches. Note
that this function potentially created new nodes. Arrival at the destination is
indicated by returning the same node that was passed as the second argument.

    updater!(tree, node, data)

Update the destination node `node`. Typically, `data` is added in some sense
to the data residing at the desitination node.
"""
function update!(tree, node, data, router!, updater!)
    while true
        dest_node = router!(tree, node)
        dest_node == node && break
        node = dest_node
    end
    updater!(tree, node, data)
    return node
end
