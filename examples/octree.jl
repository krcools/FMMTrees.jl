using FMMTrees
using StaticArrays

tree = FMMTrees.Octrees.Octree(Int[])

smallest_box_size = 0.1
root_center = SVector{3,Float64}(0,0,0)
root_size = 1.0
for i in 1:100
    target_point = rand(SVector{3,Float64})
    router! = FMMTrees.Octrees.Router(smallest_box_size, target_point)
    state = (root(tree), root_center, root_size)
    update!(tree, state, i, router!, FMMTrees.Octrees.updater!)
end

FMMTrees.print_tree(tree)

num_printed = 0
for node in FMMTrees.DepthFirstIterator(tree, root(tree))
    println(node, ", ", FMMTrees.data(tree,node))
    global num_printed += 1
end

@assert num_printed == length(tree.pbtree.nodes)
