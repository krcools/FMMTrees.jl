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
Base.IteratorSize(cv::ChildView) = Base.SizeUnknown()
AbstractTrees.children(tree::PointerBasedTree) = collect(ChildView(tree))

function Base.iterate(cv::ChildView, state=start(cv))
    done(cv) && return nothing
    return next(cv, state)
end
start(cv::ChildView) = (cv.tree.nodes[cv.tree.root].first_child)
done(cv::ChildView, idx) = (idx == -1)
function next(cv::ChildView, idx)
    value = PointerBasedTree(cv.tree.nodes, idx)
    idx = cv.tree.nodes[idx].next_sibling
    return value, idx
end

Base.iterate(cv::ChildView, s=start(cv)) = done(cv,s) ? nothing : next(cv,s)

FMMTrees.data(tree::PointerBasedTree) = tree.nodes[tree.root].data
AbstractTrees.printnode(io::IO, tree::PointerBasedTree) = show(io, data(tree))

end # module PointerBasedTrees
