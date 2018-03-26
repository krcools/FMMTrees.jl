module PointerBasedTrees

using AbstractTrees
using FMMTrees

struct Node{T}
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

struct ChildView{T} tree::T end
Base.iteratorsize(cv::ChildView) = Base.SizeUnknown()
AbstractTrees.children(tree::PointerBasedTree) = collect(ChildView(tree))

Base.start(cv::ChildView) = (cv.tree.nodes[cv.tree.root].first_child)
Base.done(cv::ChildView, idx) = (idx == -1)
function Base.next(cv::ChildView, idx)
    value = PointerBasedTree(cv.tree.nodes, idx)
    idx = cv.tree.nodes[idx].next_sibling
    return value, idx
end

FMMTrees.data(tree::PointerBasedTree) = tree.nodes[tree.root].data
AbstractTrees.printnode(io::IO, tree::PointerBasedTree) = showcompact(io, data(tree))

end # module PointerBasedTrees
