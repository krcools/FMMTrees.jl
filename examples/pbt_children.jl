using FMMTrees
using StaticArrays

struct Data
    values::Vector{Int}
end

Data(i=-1) = Data(Int[i])

const N = FMMTrees.PointerBasedTrees.Node{Data}
tree = FMMTrees.PointerBasedTrees.PointerBasedTree(
    N[
        N(Data(1),2,0,0,2),
            N(Data(2),2,5,1,3),
                N(Data(3),0,4,2,0),
                N(Data(4),0,0,2,0),
            N(Data(5),0,0,1,0)], 1)


for c in FMMTrees.children(tree, tree.root)
    println(c)
end

FMMTrees.print_tree(tree, tree.root)
