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

c = SVector{3,Float64}(1,1,1)
p = SVector{3,Float64}(0,3,-2)

s, c2, h = sector_center_size(p, c, 1.0)
s == 2
c2 == SVector{3,Float64}(0.5,1.5,0.5)
h = 0.5

function lastnonemptychild(tree,node)
    @assert FMMTrees.PointerBasedTrees.getnode(tree,node).first_child >= 1
    r = 0
    for chd in children(tree, node)
        if !(FMMTrees.PointerBasedTrees.getnode(tree,chd).first_child < 1)
            r = chd
        end
    end
    return r
end


function lastchild(tree,node)
    @assert FMMTrees.PointerBasedTrees.getnode(tree,node).first_child >= 1
    r = 0
    for chd in children(tree, node)
        r = chd
    end
    return r
end

struct Router{P,T}
    target_point::P
    smallest_box_size::T
end

function reached(tree, target, meta)
    sector, center, size = meta
    return size <= target.smallest_box_size
end

function directions(tree, target, meta)
    sector, center, size = meta
    return sector_center_size(target.target_point, center, size)
end

function isontherighttrack(tree, node, meta)
    sector, center, size = meta
    FMMTrees.data(tree, node).sector == sector
end

function newnode!(tree, node, prev, meta)
    sector, center, size = meta
    data = Data(sector)
    bef = prev < 1 ? 0 : FMMTrees.PointerBasedTrees.getnode(tree, prev).next_sibling
    FMMTrees.insert!(tree, data, parent=node, next=bef, prev=prev)
end

function FMMTrees.route!(tree::FMMTrees.PointerBasedTrees.PointerBasedTree, state, target)

    (parent, prev_fat_par, meta) = state
    @show state
    @assert prev_fat_par < 1 || FMMTrees.haschildren(tree, prev_fat_par)
    reached(tree, target, meta) && return state

    meta = directions(tree, target, meta)
    prev_fat_child = prev_fat_par < 1 ? 0 : lastnonemptychild(tree, prev_fat_par)
    prev_child = prev_fat_par < 1 ? 0 : lastchild(tree, prev_fat_par)
    for child in FMMTrees.children(tree, parent)
        isontherighttrack(tree, child, meta) && return (child, prev_fat_child, meta)
        FMMTrees.haschildren(tree, child) && (prev_fat_child = child)
        prev_child = child
    end

    child = newnode!(tree, parent, prev_child, meta)
    return child, prev_fat_child, meta
end

const N = FMMTrees.PointerBasedTrees.Node{Data}
tree = FMMTrees.PointerBasedTrees.PointerBasedTree(
    N[N(Data(), 0, 0, 0, 0)], 1)

function update!(f, tree::FMMTrees.PointerBasedTrees.PointerBasedTree, i::Int, point, sms::Float64)
    router! = Router(point, sms)
    prev_fat_child = 0
    root_sector, root_center, root_size = 0, SVector{3,Float64}(0,0,0), 1.0
    root_meta = root_sector, root_center, root_size
    root_state = (root(tree), prev_fat_child, root_meta)
    FMMTrees.update!(f, tree, root_state, i, router!)
end

using DelimitedFiles
Q = readdlm(joinpath(@__DIR__, "points.dlm"), Float64)
points = [SVector{3,Float64}(Q[i,:]) for i in axes(Q,1)]
num_points = length(points)

smallest_box_size = 0.2
f(tree, state, data) = push!(FMMTrees.data(tree, state[1]).values, data)

for i in 1:num_points
    update!(tree, i, points[i], smallest_box_size) do tree, node, data
        push!(FMMTrees.data(tree, node).values, data)
    end
end

ns = sum(length(data(tree,nd).values) for nd in FMMTrees.DepthFirstIterator(tree, tree.root))
@assert ns == num_points

FMMTrees.print_tree(tree)

num_printed = 0
for node in FMMTrees.DepthFirstIterator(tree, tree.root)
    global num_printed += 1
end
@assert num_printed == length(tree.nodes)
