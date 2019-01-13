using FMMTrees
using StaticArrays

struct Data
    sector::Int
    values::Vector{Int}
end

Data() = Data(0, Int[])

function sector_center_size(pt, ct, hs)
    hs = hs / 2
    bl = SVector([x > y for (x,y) in zip(pt,ct)]...)
    ct = SVector([b ? y+hs : y-hs for (b,y) in zip(bl,ct)]...)
    sc = sum(b ? 2^(i-1) : 0 for (i,b) in enumerate(bl))
    return sc, ct, hs
end

c = SVector{3,Float64}(1,1,1)
p = SVector{3,Float64}(0,3,-2)

s, c2, h = sector_center_size(p, c, 1.0)
s == 2
c2 == SVector{3,Float64}(0.5,1.5,0.5)
h = 0.5

function router!(tree, state)
    node, center, size, idx = state
    size <= smallest_box_size && return state
    child_sector, child_center, child_size = sector_center_size(point, center, size)
    for child in FMMTrees.children(tree, node)
        child.data.sector == child_sector && return (child, child_center, child_size)
    end
    new_child = FMMTrees.PointerBasedTrees.Node(
        Data(child_sector, Int[]), 0, node.first_child, idx, 0)
    push!(tree.nodes, new_child)
    tree.nodes[idx].first_child = length(tree.nodes)
    return new_child, child_center, child_size, length(tree.nodes)
end

function updater!(tree, state, data)
    push!(state[1].data.values, data)
end


const N = FMMTrees.PointerBasedTrees.Node{Data}
tree = FMMTrees.PointerBasedTrees.PointerBasedTree(
    N[N(Data(), 0, 0, 0, 0)], 1)

smallest_box_size = 0.1
point = SVector{3,Float64}(0,3,-2)
state = (tree.nodes[tree.root], SVector{3,Float64}(1,1,1), 1.0, tree.root)

FMMTrees.update!(tree, state, 123, router!, updater!)

FMMTrees.print_tree(tree)
