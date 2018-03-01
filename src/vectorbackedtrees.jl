struct TreeNode{T}
    num_children::Int
    data::T
end

struct VectorBackedTree{V <: (AbstractVector{N} where {N<:TreeNode})}
    nodes::V
end

Base.start(itr::ChildView{T} where {T<:VectorBackedTree}) = (0, 2) # progess, relative_index
Base.done(itr::ChildView{T} where {T<:VectorBackedTree}, state) = (state[1] == itr.tree[1].num_children)
function Base.next(itr::ChildView{T} where {T<:VectorBackedTree}, state)
    child = itr.tree[state[2]]
    newstate = (state[1] + child.num_children + 1, state[2] + child.num_children + 1)
    return view(itr.tree, state[2]:endof(itr.tree)) , newstate
end
