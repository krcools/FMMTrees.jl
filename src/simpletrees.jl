struct TreeNode{T}
    num_children::Int
    data::T
end

struct SimpleTree{V <: (AbstractVector{N} where {N<:TreeNode})}
    nodes::V
end

import ..FMMTrees: root
root(tree::SimpleTree) = tree.nodes

Base.start(itr::ChildView{T} where {T<:SimpleTree}) = (0, 2) # progess, relative_index
Base.done(itr::ChildView{T} where {T<:SimpleTree}, state) = (state[1] == itr.tree[1].num_children)
function Base.next(itr::ChildView{T} where {T<:SimpleTree}, state)
    child = itr.tree[state[2]]
    newstate = (state[1] + child.num_children + 1, state[2] + child.num_children + 1)
    return view(itr.tree, state[2]:endof(itr.tree)) , newstate
end
