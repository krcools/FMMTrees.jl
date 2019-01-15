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
    f(tree, node, level) || return
    for c in FMMTrees.children(tree, node)
        depthfirst(f, tree, c, level+1)
    end
end

function print_tree(tree, node=root(tree); maxdepth=0)
    depthfirst(tree) do tree, node, level
        print("-"^(level-1))
        print(data(tree, node))
        println()
        return level == maxdepth ? false : true
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

struct DepthFirstIterator{T,N}
    tree::T
    node::N
end

function Base.iterate(itr::DepthFirstIterator)
    chitr = children(itr.tree,itr.node)
    stack = Any[(chitr, iterate(chitr))]
    iterate(itr, stack)
end

function Base.iterate(itr::DepthFirstIterator, stack)
    isempty(stack) && return nothing
    while true
        chditr, next = last(stack)
        if next != nothing
            node, state = next
            chitr = children(itr.tree, node)
            push!(stack, (chitr, iterate(chitr)))
        else
            pop!(stack)
            isempty(stack) && return (itr.node, stack)
            chitr, (node, state) = last(stack)
            stack[end] = (chitr, iterate(chitr, state))
            return node, stack
        end
    end
end


depthfirst(tree, node=root(tree)) = DepthFirstIterator(tree,node)
