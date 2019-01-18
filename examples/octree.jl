using FMMTrees
using StaticArrays

const P = SVector{3,Float64}

tree = FMMTrees.Octrees.Octree(Int[])
points = [rand(P) for i in 1:100]

smallest_box_size = 0.1
root_center = P(0,0,0)
root_size = 1.0
for i in 1:100
    target_point = points[i]
    router! = FMMTrees.Octrees.Router(smallest_box_size, target_point)
    state = (root(tree), root_center, root_size, 1)
    update!(tree, state, i, router!, FMMTrees.Octrees.updater!)
end

FMMTrees.print_tree(tree)

num_bins = 4
bins = [P[] for i in 1:num_bins]

num_printed = 0
num_nodes = length(tree.pbtree.nodes)
for (i,node) in enumerate(FMMTrees.DepthFirstIterator(tree, root(tree)))
    println(node, ", ", FMMTrees.data(tree,node))
    b = (div((i-1) * num_bins, num_nodes) + 1
    b = i / num_nodes *
    append!(bins[b], points[FMMTrees.data(tree,node).values])
    # if i <= div(length(tree.pbtree.nodes),2)
    #     append!(first_half, points[FMMTrees.data(tree,node).values])
    # else
    #     append!(second_half, points[FMMTrees.data(tree,node).values])
    # end
    global num_printed += 1
end

@assert num_printed == length(tree.pbtree.nodes)

plot()
for i in 1:num_bins
    plot!(bins[i])
end
