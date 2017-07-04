using FMMTrees
using Base.Test

p = rand(20)
q, b = FMMTrees.clustertree(p)

children = FMMTrees.children

chdit = children(b)
state = start(children(b))


print_tree(b)
