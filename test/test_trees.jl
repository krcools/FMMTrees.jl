using FMMTrees
using Base.Test

p = rand(20)
q, b = FMMTrees.clustertree(p)

tree = FMMTrees.VectorBackedTree(b)
children = FMMTrees.children

chdit = children(tree)
state = start(children(tree))
