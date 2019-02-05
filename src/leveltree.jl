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

FMMTrees.PointerBasedTrees.data(n::HNode) = data(n.node)

struct LevelledTree{D} <: FMMTrees.PointerBasedTrees.APBTree
    nodes::Vector{HNode{D}}
    root::Int
end

FMMTrees.root(tree::LevelledTree) = tree.root

FMMTrees.PointerBasedTrees.getnode(tree::LevelledTree, idx) = tree.nodes[idx]
FMMTrees.PointerBasedTrees.nextsibling(tree::LevelledTree, node_idx) = getnode(tree, node_idx).node.next_sibling
FMMTrees.PointerBasedTrees.parent(tree::LevelledTree, node_idx) = getnode(tree, node_idx).node.parent
FMMTrees.PointerBasedTrees.firstchild(tree::LevelledTree, node_idx) = getnode(tree, node_idx).node.first_child

FMMTrees.PointerBasedTrees.setfirstchild!(node::HNode, child) = HNode(FMMTrees.PointerBasedTrees.Node(node.data, node.num_children, node.next_sibling, node.parent, child), node.height)
FMMTrees.PointerBasedTrees.setfirstchild!(tree::LevelledTree, node, child) = tree.nodes[node] = setfirstchild!(getnode(tree, node), child)

FMMTrees.PointerBasedTrees.setnextsibling!(node::HNode, next) = HNode(FMMTrees.PointerBasedTrees.Node(node.data, node.num_children, next, node.parent, node.first_child), height)
FMMTrees.PointerBasedTrees.setnextsibling!(tree::LevelledTree, node, next) = tree.nodes[node] = setnextsibling!(getnode(tree, node), next)


function FMMTrees.insert!(tree::LevelledTree, data; parent, next, prev)
    push!(tree.nodes, HNode(Node(data, 0, before, parent, 0), 0))
    fs = firstchild(tree, parent)
    if fs < 1 || fs == before
        setfirstchild!(tree, parent, length(tree.nodes))
    end
    if !(prev < 1)
        setnextsibling!(tree, prev, length(tree.nodes))
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

function route!(tree::LevelledTree, state, router)

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
