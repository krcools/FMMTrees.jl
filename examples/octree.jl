using FMMTrees
using StaticArrays
using CompScienceMeshes

const P = SVector{3,Float64}

tree = FMMTrees.Octrees.Octree(Int[])
mesh = meshsphere(1.0, 0.15)
# points = [rand(P) for i in 1:800]
points = vertices(mesh)

smallest_box_size = 0.1
root_center = P(0,0,0)
root_size = 1.0
for i in 1:length(points)
    target_point = points[i]
    router! = FMMTrees.Octrees.Router(smallest_box_size, target_point)
    state = (root(tree), root_center, root_size, 1)
    update!(tree, state, i, router!, FMMTrees.Octrees.updater!)
end

FMMTrees.print_tree(tree)

num_bins = 8
bins = [P[] for i in 1:num_bins]

num_printed = 0
num_nodes = length(tree.pbtree.nodes)
for (i,node) in enumerate(FMMTrees.DepthFirstIterator(tree, root(tree)))
    println(node, ", ", FMMTrees.data(tree,node))
    b = div((i-1)*num_bins, num_nodes) + 1
    append!(bins[b], points[FMMTrees.data(tree,node).values])
    global num_printed += 1
end

ordered_points = reduce(append!, bins)
ordered_points = [p[i] for p in ordered_points, i = 1:3]
@assert num_printed == length(tree.pbtree.nodes)

# using Plots
# plot()
# for i in 1:length(bins)
#     scatter!(bins[i])
# end
# plot!(legend=false)
