using FMMTrees
using StaticArrays
using CompScienceMeshes

const P = SVector{3,Float64}

rootnode = FMMTrees.LevelledTrees.HNode(FMMTrees.PointerBasedTrees.Node(FMMTrees.LevelledTrees.Data(0,Int[]), 0, 0, 0, 0), 0)
tree = FMMTrees.LevelledTrees.LevelledTree([rootnode], 1)
mesh = meshsphere(1.0, 0.15)
# points = [rand(P) for i in 1:800]
points = vertices(mesh)

struct OTRouter
    smallest_box_size::Float64
    target_point::P
end

FMMTrees.route!(t::FMMTrees.Octrees.Octree, state, router) = FMMTrees.Octrees.route!(t, state, router)

smallest_box_size = 0.1
root_sector = 0
root_center = P(0,0,0)
root_size = 1.0
for i in 1:length(points)
    router = OTRouter(smallest_box_size, points[i])
    root_state = root(tree), root_center, root_size, 1
    update!(tree, root_state, i, router, FMMTrees.Octrees.updater!)
end

FMMTrees.print_tree(tree)

num_bins = 8
bins = [P[] for i in 1:num_bins]

num_printed = 0
num_points = 0
num_nodes = length(tree.pbtree.nodes)
for (i,node) in enumerate(FMMTrees.DepthFirstIterator(tree, root(tree)))
    println(node, ", ", FMMTrees.data(tree,node))
    b = div((i-1)*num_bins, num_nodes) + 1
    append!(bins[b], points[FMMTrees.data(tree,node).values])
    global num_printed += 1
    global num_points += length(data(tree,node).values)
end

ordered_points = reduce(append!, bins)
ordered_points = [p[i] for p in ordered_points, i = 1:3]
@assert num_printed == length(tree.pbtree.nodes)
@show num_points
@show length(points)
@assert num_points == length(points)

# using Plots
# plot()
# for i in 1:length(bins)
#     scatter!(bins[i])
# end
# plot!(legend=false)
