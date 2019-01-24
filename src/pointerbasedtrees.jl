module PointerBasedTrees
using FMMTrees

mutable struct Node{T}
    data::T
    num_children::Int
    next_sibling::Int
    parent::Int
    first_child::Int
end

abstract type AbstractPBTree end

struct PointerBasedTree{N<:Node} <: AbstractPBTree
    nodes::Vector{N}
    root::Int
end

FMMTrees.root(tree::PointerBasedTree) = tree.root
getnode(tree::PointerBasedTree, node_idx) = tree.nodes[node_idx]

struct ChildView{T,N}
    tree::T
    node_idx::N
end

function Base.iterate(itr::ChildView)
    node = getnode(itr.tree, itr.node_idx)
    iterate(itr, node.first_child)
end

function Base.iterate(itr::ChildView, state)
    # Check if state is a valid node pointer
    state < 1 && return nothing
    # If yes, check that the parent of the corresponding node equals
    # the original parent. If not, we are looking at a distant *cousin*
    sibling_par = getnode(itr.tree, state).parent
    sibling_par != itr.node_idx && return nothing
    sibling_idx = getnode(itr.tree, state).next_sibling
    return (state, sibling_idx)
end

Base.IteratorSize(cv::ChildView) = Base.SizeUnknown()
FMMTrees.children(tree::PointerBasedTree, node=FMMTrees.root(tree)::Node) = ChildView(tree, node)
FMMTrees.data(tree::PointerBasedTree, node=FMMTrees.root(tree)) = getnode(tree, node).data
FMMTrees.haschildren(tree::PointerBasedTree, node) = (getnode(tree,node).first_child >= 1)

"""
    insert!(tree, parent, data)

Insert a node carrying `data` as the first child of `parent`
"""
function FMMTrees.insert!(tree::PointerBasedTree, parent_idx, data)
    bef = getnode(tree, parent_idx).first_child
    FMMTrees.insert!(tree, data, parent=parent_idx, before=bef, prev=0)
end


"""
    insert!(tree, data, before=bef, parent=par, prev=prev)

Insert a new node carrying data such that the next sibling of this new node
will be `bef`. The parent node and the previous node to `bef` are required
to update all links in the data structure. If `bef` is the first child of `par`,
the previous node `prev` should be pasesed as `0`. Similarly, if the new node
should be inserted as the last child, `bef` should be passed as `0`. In the
special case of a parent node without any existing children, both `bef` and
`prev` should be set to `0`.

It may seem redundant to pass also `prev`, but this is a necessity of the single
linked list implementation of `PointerBasedTree`.
"""
function FMMTrees.insert!(tree::PointerBasedTree, data; parent, before, prev)
    push!(tree.nodes, Node(data, 0, before, parent, 0))
    @assert !(parent < 1)

    fs = getnode(tree, parent).first_child
    if fs < 1 || fs == before
        getnode(tree, parent).first_child = length(tree.nodes)
    end
    if !(prev < 1)
        getnode(tree, prev).next_sibling = length(tree.nodes)
    end
    return length(tree.nodes)
end

end
