module PointerBasedTrees

# using AbstractTrees
using FMMTrees

mutable struct Node{T}
    data::T
    num_children::Int
    next_sibling::Int
    parent::Int
    first_child::Int
end

struct PointerBasedTree{N<:Node}
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
    # node.first_child < 1 && return iterate(itr, nothing)
    iterate(itr, node.first_child)
end

function Base.iterate(itr::ChildView, state)
    state < 1 && return nothing
    sibling_idx = getnode(itr.tree, state).next_sibling
    return (state, sibling_idx)
end

Base.IteratorSize(cv::ChildView) = Base.SizeUnknown()
FMMTrees.children(tree::PointerBasedTree, node=FMMTrees.root(tree)::Node) = ChildView(tree, node)

# start(cv::ChildView) = (cv.tree.nodes[cv.tree.root].first_child)
# done(cv::ChildView, idx) = (idx == -1)
# function next(cv::ChildView, idx)
#     value = PointerBasedTree(cv.tree.nodes, idx)
#     idx = cv.tree.nodes[idx].next_sibling
#     return value, idx
# end
#
# Base.iterate(cv::ChildView, s=start(cv)) = done(cv,s) ? nothing : next(cv,s)

FMMTrees.data(tree::PointerBasedTree, node=FMMTrees.root(tree)) = getnode(tree, node).data
# AbstractTrees.printnode(io::IO, tree::PointerBasedTree) = show(io, data(tree))


function insert_child!(tree::PointerBasedTree, parent_idx, data)
    parn = tree.nodes[parent_idx]
    node = FMMTrees.PointerBasedTrees.Node(data, 0, parn.first_child, parent_idx, 0)
    push!(tree.nodes, node)
    parn.first_child = length(tree.nodes)
end


end # module PointerBasedTrees
