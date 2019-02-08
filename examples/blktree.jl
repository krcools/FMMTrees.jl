using FMMTrees
using StaticArrays
const P = SVector{3,Float64}

using DelimitedFiles
Q = readdlm(joinpath(@__DIR__,"points.dlm"), Float64)
points = [SVector{3,Float64}(Q[i,:]) for i in axes(Q,1)]

root_center = P(0.5,0.5,0.5)
root_size = 0.5
root_node = FMMTrees.LevelledTrees.HNode(FMMTrees.PointerBasedTrees.Node(FMMTrees.LevelledTrees.Data(0,Int[]), 0, 0, 0, 0), 0)
tree = FMMTrees.LevelledTrees.LevelledTree([root_node], 1, root_center, root_size, Int[1])

smallest_box_size = 0.1
root_sector = 0
root_sfc_state = 1
root_depth = 1
for i in 1:length(points)
    router = FMMTrees.LevelledTrees.Router(smallest_box_size, points[i])
    root_state = root(tree), root_center, root_size, root_sfc_state, root_depth
    FMMTrees.update!(tree, root_state, i, router) do tree, node, data
        push!(FMMTrees.data(tree,node).values, data)
    end
end


bt = FMMTrees.BlockTrees.BlockTree(tree,tree)
FMMTrees.print_tree(bt)

num_levels = length(tree.levels)
num_nodes_per_level = [count(x->true,FMMTrees.LevelledTrees.LevelIterator(tree,i)) for i in 1:num_levels]

num_nodes_in_block_tree = count(x->true,FMMTrees.DepthFirstIterator(bt, FMMTrees.root(bt)))
num_leaves_in_block_tree = count(x->true,FMMTrees.leaves(bt, FMMTrees.root(bt)))

@assert num_nodes_in_block_tree == sum(m^2 for m in num_nodes_per_level)
@assert num_leaves_in_block_tree == num_nodes_per_level[end]^2
