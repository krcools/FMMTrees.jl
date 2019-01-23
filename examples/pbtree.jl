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

function router!(tree::FMMTrees.PointerBasedTrees.PointerBasedTree, state)
    parent, center, size, prev_par = state
    @assert prev_par < 1 || FMMTrees.haschildren(tree, prev_par)

    # reached(tree, state, captured) && return state
    size <= smallest_box_size && return state

    child_sector, child_center, child_size = sector_center_size(target_point, center, size)

    prev_fat_child = prev_par < 1 ? 0 : lastnonemptychild(tree, prev_par)
    prev_child = prev_par < 1 ? 0 : lastchild(tree, prev_par)
    for child in FMMTrees.children(tree, parent)
        # belongs(tree, child, state, meta) && return newstate(tree, child, state, meta, prev_fat_child)
        d = FMMTrees.data(tree,child)
        if d.sector == child_sector
            return (child, child_center, child_size, prev_fat_child)
        end
        FMMTrees.haschildren(tree, child) && (prev_fat_child = child)
        prev_child = child
    end
    # node = newnode!(tree, state, meta, prev_child)
    data = Data(child_sector, Int[])
    bef = prev_child < 1 ? 0 : FMMTrees.PointerBasedTrees.getnode(tree, prev_child).next_sibling
    child = FMMTrees.insert!(tree, data, parent=parent, before=bef, prev=prev_child)
    # return newstate(tree, child, state, meta, prev_fat_child)
    return child, child_center, child_size, prev_fat_child
end

function updater!(tree, state, data)
    push!(tree.nodes[state[1]].data.values, data)
end


const N = FMMTrees.PointerBasedTrees.Node{Data}
tree = FMMTrees.PointerBasedTrees.PointerBasedTree(
    N[N(Data(), 0, 0, 0, 0)], 1)

smallest_box_size = 0.2
root_center = SVector{3,Float64}(0,0,0)
root_size = 1.0
root_prev = 0
using DelimitedFiles
Q = readdlm("points.dlm", Float64)
points = [SVector{3,Float64}(Q[i,:]) for i in axes(Q,1)]
points = [ rand(SVector{3,Float64}) for i in 1:100 ]
num_points = length(points)
for i in 1:num_points
    global target_point = points[i] #rand(SVector{3,Float64})
    @show target_point
    state = (tree.root, root_center, root_size, root_prev)
    FMMTrees.update!(tree, state, i, router!, updater!)
end

ns = sum(length(data(tree,nd).values) for nd in FMMTrees.DepthFirstIterator(tree, tree.root))
@assert ns == num_points


FMMTrees.print_tree(tree)

num_printed = 0
for node in FMMTrees.DepthFirstIterator(tree, tree.root)
    println(node, ", ", FMMTrees.data(tree,node))
    global num_printed += 1
end

@assert num_printed == length(tree.nodes)
