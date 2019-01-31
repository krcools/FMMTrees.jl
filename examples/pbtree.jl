using FMMTrees
using StaticArrays

struct Data
    sector::Int
    values::Vector{Int}
end

Data(s=0) = Data(s, Int[])

function sector_center_size(pt, ct, hs)
    hs = hs / 2
    bl = SVector([x > y for (x,y) in zip(pt,ct)]...)
    ct = SVector([b ? y+hs : y-hs for (b,y) in zip(bl,ct)]...)
    sc = sum(b ? 2^(i-1) : 0 for (i,b) in enumerate(bl))
    return sc, ct, hs
end

struct Router{P,T}
    target_point::P
    smallest_box_size::T
end

function FMMTrees.PointerBasedTrees.reached(tree, target, meta)
    sector, center, size = meta
    return size <= target.smallest_box_size
end

function FMMTrees.PointerBasedTrees.directions(tree, target, meta)
    sector, center, size, sfc_state = meta
    sector, center, size = sector_center_size(target.target_point, center, size)
    sfc_pos = FMMTrees.Octrees.hilbert_positions[sfc_state][sector+1] + 1
    sfc_state = FMMTrees.Octrees.hilbert_states[sfc_state][sector+1] + 1
    return sector, center, size, sfc_state, sfc_pos
end

function FMMTrees.PointerBasedTrees.isontherighttrack(tree, node, meta)
    sector, center, size = meta
    FMMTrees.data(tree, node).sector == sector
end

function FMMTrees.PointerBasedTrees.reachedposamongsiblings(tree, child, par_meta)
    _1, _2, _3, par_sfc_state = par_meta
    child_sector = FMMTrees.data(tree, child).sector
    child_pos = FMMTrees.Octrees.hilbert_positions[par_sfc_state][child_sector+1] + 1
    target_pos = FMMTrees.Octrees.hilbert_positions[par_sfc_state][target_sector+1] + 1
    return target_pos < child_pos
end

function FMMTrees.PointerBasedTrees.newnode!(tree, node, prev, bef, meta)
    sector, center, size, sfc_state = meta
    data = Data(sector)
    # bef = prev < 1 ? 0 : FMMTrees.PointerBasedTrees.getnode(tree, prev).next_sibling
    @show node, bef, prev
    FMMTrees.insert!(tree, data, parent=node, before=bef, prev=prev)
end

function updater!(tree, state, data)
    node_idx = state[1]
    push!(tree.nodes[node_idx].data.values, data)
end

const N = FMMTrees.PointerBasedTrees.Node{Data}
tree = FMMTrees.PointerBasedTrees.PointerBasedTree(
    N[N(Data(), 0, 0, 0, 0)], 1)

function FMMTrees.update!(f, tree::FMMTrees.PointerBasedTrees.PointerBasedTree, i, point, sms)
    router! = Router(point, sms)
    prev_fat_child = 0
    root_sector, root_center, root_size = 0, SVector{3,Float64}(0.5,0.5,0.5), 0.5
    root_sfc_state = 1
    root_sfc_pos = 1
    root_meta = root_sector, root_center, root_size
    root_state = (root(tree), prev_fat_child, root_sfc_state, root_meta)
    FMMTrees.update!(tree, root_state, i, router!, f)
end

using DelimitedFiles
Q = readdlm("points.dlm", Float64)
points = [SVector{3,Float64}(Q[i,:]) for i in axes(Q,1)]
# points = [ rand(SVector{3,Float64}) for i in 1:10 ]
num_points = length(points)

smallest_box_size = 0.2
for i in 1:num_points
    FMMTrees.update!(tree, i, points[i], smallest_box_size) do tree, state, data
        push!(FMMTrees.data(tree, state[1]).values, data)
    end
end

error("stop")

ns = sum(length(data(tree,nd).values) for nd in FMMTrees.DepthFirstIterator(tree, tree.root))
@assert ns == num_points

FMMTrees.print_tree(tree)

num_printed = 0
for node in FMMTrees.DepthFirstIterator(tree, tree.root)
    global num_printed += 1
end
@assert num_printed == length(tree.nodes)
