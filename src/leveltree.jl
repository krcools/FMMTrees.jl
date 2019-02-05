import FMMTrees.PointerBasedTrees
import FMMTrees.PointerBasedTrees: PointerBasedTree

struct HNode{D}
    node::PointerBasedTrees.Node{D}
    height::Int
end

PointerBasedTrees.data(n::HNode) = data(n.node)

struct LevelledTree{D} <: PointerBasedTrees.APBTree
    nodes::Vector{HNode{D}}
    root::Int
end

FMMTrees.root(tree::LevelledTree) = tree.root

PointerBasedTrees.getnode(tree::LevelledTree, idx) = tree.nodes[idx]
PointerBasedTrees.nextsibling(tree::LevelledTree, node_idx) = getnode(tree, node_idx).node.next_sibling
PointerBasedTrees.parent(tree::LevelledTree, node_idx) = getnode(tree, node_idx).node.parent
PointerBasedTrees.firstchild(tree::LevelledTree, node_idx) = getnode(tree, node_idx).node.first_child

PointerBasedTrees.setfirstchild!(node::HNode, child) = HNode(PointerBasedTrees.Node(node.data, node.num_children, node.next_sibling, node.parent, child), node.height)
PointerBasedTrees.setfirstchild!(tree::LevelledTree, node, child) = tree.nodes[node] = setfirstchild!(getnode(tree, node), child)

PointerBasedTrees.setnextsibling!(node::HNode, next) = HNode(PointerBasedTrees.Node(node.data, node.num_children, next, node.parent, node.first_child), height)
PointerBasedTrees.setnextsibling!(tree::LevelledTree, node, next) = tree.nodes[node] = setnextsibling!(getnode(tree, node), next)


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
