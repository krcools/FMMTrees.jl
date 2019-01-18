module Octrees

using FMMTrees
using StaticArrays

struct Data{T}
    sector::Int
    values::Vector{T}
end

struct Octree{T}
    pbtree::FMMTrees.PointerBasedTrees.PointerBasedTree{FMMTrees.PointerBasedTrees.Node{Data{T}}}
    function Octree(value::Vector{T}) where {T}
        N = FMMTrees.PointerBasedTrees.Node{Data{T}}
        new{T}(
            FMMTrees.PointerBasedTrees.PointerBasedTree(
                N[N(Data(0,value),0,0,0,0)], 1))
    end
end

FMMTrees.data(tree::Octree, node) = data(tree.pbtree, node)
FMMTrees.children(tree::Octree, node) = children(tree.pbtree, node)
FMMTrees.insert!(tree::Octree, node, data) = FMMTrees.insert!(tree.pbtree, node, data)
FMMTrees.insert!(tree::Octree, data; before, parent, prev) = FMMTrees.insert!(tree.pbtree, data, before=before, parent=parent, prev=prev)
FMMTrees.root(tree::Octree) = root(tree.pbtree)

const hilbert_states = [
    [1, 2, 3, 2, 4, 5, 3, 5],
    [2, 6, 0, 7, 8, 8, 0, 7],
    [0, 9,10, 9, 1, 1,11,11],
    [6, 0, 6,11, 9, 0, 9, 8],
    [11,11, 0, 7, 5, 9, 0, 7],
    [4, 4, 8, 8, 0, 6,10, 6],
    [5, 7, 5, 3, 1, 1,11,11],
    [6, 1, 6,10, 9, 4, 9,10],
    [10, 3, 1, 1,10, 3, 5, 9],
    [4, 4, 8, 8, 2, 7, 2, 3],
    [7, 2,11, 2, 7, 5, 8, 5],
    [10, 3, 2, 6,10, 3, 4, 4]]

const hilbert_positions = [
    [0,1,3,2,7,6,4,5],
    [0,7,1,6,3,4,2,5],
    [0,3,7,4,1,2,6,5],
    [2,3,1,0,5,4,6,7],
    [4,3,5,2,7,0,6,1],
    [6,5,1,2,7,4,0,3],
    [4,7,3,0,5,6,2,1],
    [6,7,5,4,1,0,2,3],
    [2,5,3,4,1,6,0,7],
    [2,1,5,6,3,0,4,7],
    [4,5,7,6,3,2,0,1],
    [6,1,7,0,5,2,4,3]]

function sector_center_size(pt, ct, hs)
    hs = hs / 2
    bl = SVector([x > y for (x,y) in zip(pt,ct)]...)
    ct = SVector([b ? y+hs : y-hs for (b,y) in zip(bl,ct)]...)
    sc = sum(b ? 2^(i-1) : 0 for (i,b) in enumerate(bl))
    return sc, ct, hs
end


struct Router{T,P}
    smallest_box_size::T
    target_point::P
end

function (f::Router)(tree, state)

    point = f.target_point
    smallest_box_size = f.smallest_box_size

    node_idx, center, size, sfc_state = state
    size <= smallest_box_size && return state
    child_sector, child_center, child_size = sector_center_size(point, center, size)
    child_sfc_state = hilbert_states[sfc_state][child_sector+1] + 1
    for child in FMMTrees.children(tree, node_idx)
        d = FMMTrees.data(tree,child)
        if d.sector == child_sector
            return (child, child_center, child_size, child_sfc_state)
        end
    end
    data = Data(child_sector, Int[])
    new_node_idx = FMMTrees.insert!(tree, (node_idx, sfc_state), data)
    return new_node_idx, child_center, child_size, child_sfc_state
end


function updater!(tree, (node, center, size), data)
    push!(FMMTrees.data(tree, node).values, data)
end


function FMMTrees.insert!(tree::Octree, (parent, sfc_state), data)

    target_pos = hilbert_positions[sfc_state][data.sector+1] + 1
    prev = 0
    found = false
    for child in children(tree, parent)
        child_sector = FMMTrees.data(tree, child)
        child_pos = hilbert_positions[sfc_state][child_sector+1] + 1
        if child_pos > target_pos
            return FMMTrees.insert!(tree, data, before=child, prev=prev, parent=parent)
        end
        prev = child
    end

    return FMMTrees.insert!(tree, data, before=0, prev=prev, parent=parent)
end


end # module
