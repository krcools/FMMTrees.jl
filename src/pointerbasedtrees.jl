module PointerBasedTrees

using ..FMMTrees

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

import FMMTreees: root, children, insertchild
root(tree::PointerBasedTree) = tree.nodes[root]
function insertchild(root, node)
    
end

end # module PointerBasedTrees
