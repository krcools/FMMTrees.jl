module  LevelledTrees

import FMMTrees

struct Data{T}
    sector::Int
    values::Vector{T}
end

struct HNode{D}
    node::FMMTrees.PointerBasedTrees.Node{D}
    height::Int
end

# FMMTrees.PointerBasedTrees.data(n::HNode) = data(n.node)

struct LevelledTree{D} <: FMMTrees.PointerBasedTrees.APBTree
    nodes::Vector{HNode{D}}
    root::Int
end

FMMTrees.root(tree::LevelledTree) = tree.root
FMMTrees.data(tree::LevelledTree, node=FMMTrees.root(tree)) = tree.nodes[node].node.data

FMMTrees.PointerBasedTrees.nextsibling(tree::LevelledTree, node_idx) = tree.nodes[node_idx].node.next_sibling
FMMTrees.PointerBasedTrees.parent(tree::LevelledTree, node_idx) = tree.nodes[node_idx].node.parent
FMMTrees.PointerBasedTrees.firstchild(tree::LevelledTree, node_idx) = tree.nodes[node_idx].node.first_child

_setfirstchild!(node::HNode, child) = HNode(FMMTrees.PointerBasedTrees.Node(node.node.data, node.node.num_children, node.node.next_sibling, node.node.parent, child), node.height)
FMMTrees.PointerBasedTrees.setfirstchild!(tree::LevelledTree, node, child) = tree.nodes[node] = _setfirstchild!(tree.nodes[node], child)

_setnextsibling!(node::HNode, next) = HNode(FMMTrees.PointerBasedTrees.Node(node.node.data, node.node.num_children, next, node.node.parent, node.node.first_child), node.height)
FMMTrees.PointerBasedTrees.setnextsibling!(tree::LevelledTree, node, next) = tree.nodes[node] = _setnextsibling!(tree.nodes[node], next)


function FMMTrees.insert!(tree::LevelledTree, data; parent, next, prev)
    push!(tree.nodes, HNode(FMMTrees.PointerBasedTrees.Node(data, 0, next, parent, 0), 0))
    fs = FMMTrees.PointerBasedTrees.firstchild(tree, parent)
    if fs < 1 || fs == next
        FMMTrees.PointerBasedTrees.setfirstchild!(tree, parent, length(tree.nodes))
    end
    if !(prev < 1)
        FMMTrees.PointerBasedTrees.setnextsibling!(tree, prev, length(tree.nodes))
    end
    return length(tree.nodes)
end

function sector_center_size(pt, ct, hs)
    hs = hs / 2
    bl = pt .> ct
    ct = ifelse.(bl, ct.+hs, ct.-hs)
    sc = sum(b ? 2^(i-1) : 0 for (i,b) in enumerate(bl))
    return sc, ct, hs
end

struct Router{T,P}
    smallest_box_size::T
    target_point::P
end

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

function FMMTrees.route!(tree::LevelledTree, state, router)

    point = router.target_point
    smallest_box_size = router.smallest_box_size

    node_idx, center, size, sfc_state = state
    size <= smallest_box_size && return state
    target_sector, target_center, target_size = sector_center_size(point, center, size)
    target_pos = hilbert_positions[sfc_state][target_sector+1] + 1
    target_sfc_state = hilbert_states[sfc_state][target_sector+1] + 1
    prev_child, next_child = 0, 0
    for child in FMMTrees.children(tree, node_idx)
        child_sector = FMMTrees.data(tree,child).sector
        child_pos = hilbert_positions[sfc_state][child_sector+1]+1
        target_pos < child_pos  && (next_child = child; break)
        if child_sector == target_sector
            return (child, target_center, target_size, target_sfc_state)
        end
        prev_child = child
    end
    data = Data(target_sector, Int[])
    new_node_idx = FMMTrees.insert!(tree, data, next=next_child, prev=prev_child, parent=node_idx)
    return new_node_idx, target_center, target_size, target_sfc_state
end

end # module LevelledTrees
