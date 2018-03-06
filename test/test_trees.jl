using FMMTrees
using Base.Test

# p = rand(20)
# q, b = FMMTrees.clustertree(p)

N = FMMTrees.TreeNode{Int}
nodes = N[
    N(6,-1),
        N(2,-1),
            N(0,-1),
            N(0,-1),
        N(4,-1),
            N(0,-1),
            N(2,-1),
                N(0,-1),
                N(0,-1)]

tree = FMMTrees.SimpleTree(nodes)
children = FMMTrees.children

chdit = children(tree, root(tree))
state = start(children(tree, root(tree)))
